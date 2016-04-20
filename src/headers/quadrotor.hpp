/*************************************************************************
    > File Name: quadrotor.hpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Tue 21 Jul 2015 11:10:23 PM EDT
 ************************************************************************/
#ifndef QUADROTOR
#define QUADROTOR

#include<list>
#include<iostream>
#include<cmath>
#include<string>
#include<functional>

#include<Eigen/Dense>
#include<mat.h>
#include<liegroup.hpp>
#include<quadcost.hpp>

#include<quadopt.hpp>
struct quadrotor
{
	public:
		struct System;
		struct State;
		struct DState;
	
//	protected:
		static quadopt optProb;


	public:
		struct System
		{
			const Mat3 I; // moment of inertia 
			const double m; // mass

			const double d; // displacement of the motor
			const double km; // coefficient to balance the torque of the motor
			const double kt; // coefficient to generate lift force

			const Mat3 I_inv; // inverse of moment inertia

			System(double Ix, double Iy, double Iz, double m_, double d_, double km_, double kt_);
		};

	//**************************************************************************************************************************
		// Definition of quadrotor state and dynamics 
	// *************************************************************************************************************************
		
		struct State // quadrotor state
		{
			Mat4 g; // SE(3) to represent orientation and velocity
			Vec6 v; // body twist

			State();
			State(Mat4 g0, Vec6 v0);
			State(Mat3 R0, Vec3 x0, Vec3 w0, Vec3 v0);

			State update(const DState & dstate, double h); // compute next state by Euler method
			static void  save(const std::list<State> &, std::string);
		};
		
		struct DState // time derivative of quadrotor state
		{
			Vec6 v; // body twist
			Vec6 a; // body twist accelebration;

			DState(const System & sys, const State & state, Vec4 u); // compute body velocity
		};

	public:
		typedef State Ref; // references used to compute cost function

		static constexpr double g=9.8; // accelebration of gravity

		static constexpr size_t M=12; // dimension of the configuration space
		static constexpr size_t N=4;  // dimension of control inputs

		typedef Eigen::Matrix<double, M,1> V;
		typedef Eigen::Matrix<double, N,1> U;

		static Eigen::Matrix<double,M,1> diff(State const & state, Ref const & ref)
		{
			return (Vec12()<<SE3::log(SE3::inverse(ref.g)*state.g), state.v-ref.v).finished();
		}

		static Eigen::Matrix<double,M,1> f(System const & sys, State const &state)
		{
			Eigen::Matrix<double,M,1> f;

			Mat3 R=state.g.block(0,0,3,3);
			Vec3 omega=state.v.head(3);

			double g=quadrotor::g;

			f.block(0,0,6,1)=state.v;
			f.block(6,0,3,1)=-sys.I_inv*SO3::hat(omega)*sys.I*omega;
			f.block(9,0,3,1)=-SO3::hat(omega)*state.v.tail(3)-g*R.row(2).transpose();

			return f;
		}

		static Eigen::Matrix<double,M,N> h(System const & sys, State const & state)
		{
			Eigen::Matrix<double,quadrotor::M,quadrotor::N> H=Eigen::Matrix<double,quadrotor::M,quadrotor::N>::Zero();

			H.block(6,0,3,3)=sys.I_inv.block(0,0,3,3);

			H.block(9,3,3,1)=SO3::e[2]/sys.m;

			const static Eigen::Matrix<double,4,4> M=(Eigen::Matrix<double,4,4>()<<0,1,0,-1,-1,0,1,0,1,-1,1,-1,1,1,1,1).finished();

			Eigen::Matrix<double,4,4> MF=M;

			MF.block(0,0,2,4)*=sys.kt*sys.d;
			MF.block(2,0,1,4)*=sys.km;
			MF.block(3,0,1,4)*=sys.kt;
			
			H=H*MF;

			return H;
		}

		static Eigen::Matrix<double, M, M> Dgf(System const & sys, State const & state, U const & u)
		{
			Vec3 omega=state.v.head(3);
			Mat3 R=state.g.block(0,0,3,3);

			Eigen::Matrix<double, M, M> df=Eigen::Matrix<double, M, M>::Zero();

			df.block(0,6,6,6)=Mat6::Identity();
			df.block(6,6,3,3)=sys.I_inv*(SO3::hat(sys.I*omega)-SO3::hat(omega)*sys.I);
			df.block(9,0,3,3)=-g*SO3::hat(R.row(2).transpose());
			df.block(9,6,3,3)=SO3::hat(state.v.tail(3));
			df.block(9,9,3,3)=-SO3::hat(omega);

			return df;	
		}

		static Eigen::Matrix<double,M,N> Duf(System const & sys, State const & state, U const & u)
		{
			return h(sys,state);
		}

		static Eigen::Matrix<double, M, M> ad(System const &sys, State const &state, U const & u)
		{
			Eigen::Matrix<double, M, M> ad=Eigen::Matrix<double,M,M>::Zero();
			ad.block(0,0,6,6)=SE3::ad(state.v);

			return ad;
		}

	protected:
		static void iLQG_back(           System const & sys,                       double      dt, 
				               std::list<State> const & list_state0, std::list<U> const & list_u0, std::list<Ref> const & list_ref, 
							              Mat12 const &           M,         Mat4 const &       R,          Mat12 const &       Mf, 
							               Vec4 const &        umin,         Vec4 const &    umax,         double const &     alpha, 
							   std::list<Eigen::Matrix<double,4,12> > & list_K,         std::list<Vec4> & list_k)
		{
			list_K.clear();
			list_k.clear();

			std::list<State>::const_reverse_iterator rit_state=list_state0.crbegin();
			std::list<Ref>::const_reverse_iterator rit_ref=list_ref.crbegin();
			std::list<U>::const_reverse_iterator rit_u=list_u0.crbegin();

			Vec12 dg=(Vec12()<<SE3::log(SE3::inverse(rit_ref->g)*rit_state->g),
							   rit_state->v-rit_ref->v).finished();

			Mat12 Vxx=quadcost::Lxx(Mf,dg);
			Vec12 Vx=quadcost::Lx(Mf,dg);

			rit_state++;
			rit_ref++;
			
			quadopt::Params params;
			Eigen::SelfAdjointEigenSolver<Mat4> es;

			while(rit_state!=list_state0.crend())
			{
				Vec4 const & u=*rit_u;

				State const & state=*rit_state;
				Vec6 const & v=rit_state->v;

				Vec12 dg=(Vec12()<<SE3::log(SE3::inverse(rit_ref->g)*rit_state->g),
									rit_state->v-rit_ref->v).finished();

				Vec4 Lu=quadcost::Lu(R,u)*dt;
				Mat4 Luu=quadcost::Luu(R,u)*dt;
				Vec12 Lx=quadcost::Lx(M,dg)*dt;
				Mat12 Lxx=quadcost::Lxx(M,dg)*dt;

				Mat6 Ad=SE3::Ad(SE3::exp(-v*dt));
				Mat6 dexp=SE3::dexp(v*dt);
				
				Mat12 dgF=quadrotor::Dgf(sys,state,u);
				dgF.block(0,0,6,12)=dexp*dgF.block(0,0,6,12);

				Eigen::Matrix<double,12,4> duF=quadrotor::Duf(sys,state,u);
				duF.block(0,0,6,4)=dexp*duF.block(0,0,6,4);

				Mat12 A=Mat12::Identity()+dgF*dt;
				A.block(0,0,6,12)=Ad*A.block(0,0,6,12);
				Mat12 At=A.transpose();

				Eigen::Matrix<double,12,4> B=duF*dt;	
				Eigen::Matrix<double,4,12> Bt=B.transpose();

				Vec12 Qx=Lx+At*Vx;
				Vec4 Qu=Lu+Bt*Vx;
				Mat12 Qxx=Lxx+At*Vxx*A;
				Eigen::Matrix<double,4,12> Qux=Bt*Vxx*A;
				Eigen::Matrix<double,12,4> Qxu=Qux.transpose();
				Mat4 Quu=Luu+Bt*Vxx*B;

				es.compute(Quu);
				
				Vec4 evals=es.eigenvalues();
				Mat4 V=es.eigenvectors();
				
				for(int i=0;i<4;i++)
					if(evals(i)<0)
						evals(i)=0;
					else
						break;
				
				evals+=alpha*Vec4::Ones();
				Mat4 D=evals.asDiagonal();
				Quu=V*D*V.transpose();
				Mat4 Quuinv=V*(Vec4()<<1.0/evals(0),1.0/evals(1),1.0/evals(2),1.0/evals(3)).finished().asDiagonal()*V.transpose();

				Vec4 dumax=umax-u;
				Vec4 dumin=umin-u;
				Vec4 k=-Quuinv*Qu;
				Eigen::Matrix<double,4,12> K=-Quuinv*Qux;

				if(k(0)>dumax(0) || k(1)>dumax(1) || k(2)>dumax(2) || k(3)>dumax(3) ||
				   k(0)<dumin(0) || k(1)<dumin(1) || k(2)<dumin(2) || k(3)<dumin(3))
				{
					Vec4 x=list_k.empty()? Vec4::Zero():list_k.back();

					params.xlow=dumin;
					params.xupp=dumax;

					Vec4 xmul=Vec4::Zero();
					Eigen::Matrix<int,4,1> xstate=Eigen::Matrix<int,4,1>::Zero();

					double F=(0.5*x.transpose()*Quu*x+x.transpose()*Qu)(0);
					double Fmul=0;
					int Fstate=0;
					int result=optProb.fmin(Quu,     Qu,   params,
											  x,   xmul,   xstate,
											  F,   Fmul,   Fstate);
					for(int i=0;i<4;i++)
						if(x(i)>dumax(i))
							x(i)=dumax(i);
						else
							if(x(i)<dumin(i))
								x(i)=dumin(i);

					if(result>=2)
					{
						for(int i=0;i<4;i++)
							if(k(i)>dumax(i))
								k(i)=dumax(i);
							else
								if(k(i)<dumin(i))
									k(i)=dumin(i);

						if(F<(0.5*k.transpose()*Quu*k+k.transpose()*Qu)(0))
							if(x(0)<=dumax(0) && x(1)<=dumax(1) && x(2)<=dumax(2) && x(3)<=dumax(3) &&
							   x(0)>=dumin(0) && x(1)>=dumin(1) && x(2)>=dumin(2) && x(3)>=dumin(3))
								k=x;
					}

				}
				
				for(int i=0;i<4;i++)
					if(fabs(dumax(i)-k(i))<1e-3 || fabs(k(i)-dumin(i))<1e-3)
						K.row(i)=Eigen::Matrix<double,1,12>::Zero();

				list_k.push_front(k);
				list_K.push_front(K);

				Eigen::Matrix<double,12,4> Kt=K.transpose();

				Vxx=Qxx+Kt*Quu*K+Kt*Qux+Qxu*K;	
				Vx=Qx+Kt*(Quu*k+Qu)+Qxu*k;

				std::cout<<k.transpose()<<std::endl;
				std::cout<<K<<std::endl;

				rit_u++;
				rit_state++;
				rit_ref++;
			}
		}

		static double iLQG_forward(          System const &         sys,               double      dt,                 size_t        N, 
				                   std::list<State> const & list_state0, std::list<U> const & list_u0, std::list<Ref> const & list_ref, 
								              Mat12 const &           M,         Mat4 const &       R,          Mat12 const &       Mf, 
								               Vec4 const &        umin,         Vec4 const &    umax, 
								  std::list<Eigen::Matrix<double,4,12> > const & list_K,     std::list<Vec4> const & list_k, 
								                              std::list<State> & list_state,          std::list<U> & list_u)
		{
			list_state.clear();
			list_u.clear();

			State state=list_state0.front();
			list_state.push_back(state);

			std::list<U>::const_iterator itr_u0=list_u0.begin();
			std::list<State>::const_iterator itr_state0=list_state0.begin();
			std::list<State>::const_iterator itr_ref=list_ref.begin();
			std::list<U>::const_iterator itr_k=list_k.begin();
			std::list<Eigen::Matrix<double,4,12> >::const_iterator itr_K=list_K.begin();
		
			double cost=0;

			for(int n=1; n<N; n++)
			{
				Vec12 dg=(Vec12()<<SE3::log(SE3::inverse(itr_ref->g)*state.g),
									state.v-itr_ref->v).finished();


				Vec12 error=(Vec12()<<SE3::log(SE3::inverse(itr_state0->g)*state.g),
									state.v-itr_state0->v).finished();
				
				Vec4 u=*itr_u0+*itr_k+*itr_K*error;

				cost+=quadcost::L(M,R,dg,u)*dt;
				
				for(int i=0;i<4;i++)
					if(u(i)>umax(i))
						u(i)=umax(i);
					else
						if(u(i)<umin(i))
							u(i)=umin(i);

				DState dstate(sys,state,u);
				state=state.update(dstate,dt);

				list_state.push_back(state);
				list_u.push_back(u);

				itr_k++;
				itr_K++;
				itr_ref++;
				itr_state0++;
				itr_u0++;
			}

			Vec12 dg=(Vec12()<<SE3::log(SE3::inverse(itr_ref->g)*state.g),
								state.v-itr_ref->v).finished();

			cost+=quadcost::L(Mf,Mat4::Zero(),dg,Vec4::Zero());

			return cost;
		}

	public:
		static bool iLQG(System const &    sys,         double     dt,                 double        T, 
				          State const & state0, std::list<U> & list_U, std::list<Ref> const & list_ref, 
						  Mat12 const &      M,   Mat4 const &      R,          Mat12 const &       Mf, 
						   Vec4 const &   umin,   Vec4 const &    umax, 
					     double const &rel_error=0.02,   size_t const & itr_max=10, 
						 double const & alpha=1, double const & coeff=1.6)
		{
			size_t N=size_t(T/dt+0.5);

			if(list_ref.size()<N)
			{
				std::cout<<"ERROR: Not enough references to compute cost functions."<<std::endl;
				return false;
			}


			if(list_U.size()<N-1)
			{
				std::list<U> list_Uf(N-list_U.size()-1, Vec4::Zero());
				list_U.splice(list_U.end(),list_Uf);
			}
			
			std::list<State> list_state0;
			std::list<U> list_u0;

			State state=state0;
			list_state0.push_back(state);

			std::list<U>::const_iterator itr_u=list_U.begin();
			std::list<State>::const_iterator itr_ref=list_ref.begin();

			double cost=0;

			for(int n=1;n<N;n++)
			{
				Vec12 dg=(Vec12()<<SE3::log(SE3::inverse(itr_ref->g)*state.g),
									state.v-itr_ref->v).finished();

				std::cout<<dg.transpose()<<std::endl;

				cost+=quadcost::L(M,R,dg,*itr_u)*dt;

				DState dstate(sys,state,*itr_u);
				state=state.update(dstate,dt);

				list_state0.push_back(state);
				list_u0.push_back(*itr_u);

				itr_u++;
				itr_ref++;
			}
		
			for(size_t n=0;n<itr_max;n++)
			{
				std::list<Eigen::Matrix<double,4,12> > list_K;
				std::list<Vec4> list_k;

				double al=alpha;
				iLQG_back(sys,dt,list_state0,list_u0,list_ref,M,R,Mf,umin,umax,al,list_K,list_k);

				std::list<State> list_state;
				std::list<U> list_u;
				double cost_new=iLQG_forward(        sys,      dt,        N, 
						                     list_state0, list_u0, list_ref, 
											           M,       R,       Mf, 
													umin,    umax, 
												  list_K,  list_k, 
											  list_state,  list_u);

				if(fabs(cost_new-cost)/cost<rel_error)
					break;

				std::cout<<cost<<" "<<cost_new<<std::endl;

				if(cost_new<cost)
				{
					list_u0=list_u;
					list_state0=list_state;

//					std::cout<<cost<<" "<<cost_new<<std::endl;
					
					if(al>0.01)
						al/=coeff;
				}
				else
				{
					al*=coeff;

					if(al>50)
						break;
				}
			}

			if(alpha<50)
				list_U=list_u0;
		}


};

quadrotor::System::System(double Ix, double Iy, double Iz, double m_, double d_, double km_, double kt_): I((Mat3()<<Ix,0,0,0,Iy,0,0,0,Iz).finished()), m(m_), d(d_), km(km_), kt(kt_), I_inv((Mat3()<<1.0/Ix, 0, 0, 0, 1.0/Iy,0,0, 0, 1.0/Iz).finished())
{
}


quadrotor::State::State():g(Mat4::Identity()), v(Vec6::Zero())
{
}

quadrotor::State::State(Mat4 g0, Vec6 v0):g(g0), v(v0)
{
}

quadrotor::State::State(Mat3 R0, Vec3 x0, Vec3 w0, Vec3 v0): g((Mat4()<<R0,x0,0,0,0,1).finished()), v((Vec6()<<w0,v0).finished())
{
}

quadrotor::State quadrotor::State::update(const DState & dstate, double h)
{

	const static Eigen::Matrix<double,1,4> row((Eigen::Matrix<double,1,4>()<<0,0,0,1).finished());

	State state_next;

	state_next.g=g*SE3::exp(dstate.v*h);
	state_next.v=v+dstate.a*h;

	return state_next;
}

void quadrotor::State::save(const std::list<State> & states, std::string path)
{
	MATFile *result;

	mxArray *R;
	mxArray *xq;
	mxArray *omega;
	mxArray *vq;

	void *p;
	
	mwSize *dims=new mwSize[3];
	
	result=matOpen(path.c_str(),"w");

	std::list<State>::const_iterator it_state;

	it_state=states.cbegin();
	dims[0]=3;
	dims[1]=3;
	dims[2]=states.size();

	R=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);


	for(p=mxGetPr(R);it_state!=states.cend();it_state++)
	{
		Mat3 R=it_state->g.block(0,0,3,3);
		memcpy(p,R.data(),sizeof(double)*dims[0]*dims[1]);
		p+=sizeof(double)*dims[0]*dims[1];
	}

	matPutVariable(result,"R",R);
	mxDestroyArray(R);

	it_state=states.cbegin();
	dims[0]=3;
	dims[1]=1;
	dims[2]=states.size();

	xq=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);


	for(p=mxGetPr(xq);it_state!=states.cend();it_state++)
	{
		memcpy(p,it_state->g.block(0,3,3,1).data(),sizeof(double)*dims[0]*dims[1]);
		p+=sizeof(double)*dims[0]*dims[1];
	}

	matPutVariable(result,"xq",xq);
	mxDestroyArray(xq);

	it_state=states.cbegin();
	dims[0]=3;
	dims[1]=1;
	dims[2]=states.size();

	vq=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);


	for(p=mxGetPr(vq);it_state!=states.cend();it_state++)
	{
		memcpy(p,it_state->v.tail(3).data(),sizeof(double)*dims[0]*dims[1]);
		p+=sizeof(double)*dims[0]*dims[1];
	}

	matPutVariable(result,"vq",vq);
	mxDestroyArray(vq);

	it_state=states.cbegin();
	dims[0]=3;
	dims[1]=1;
	dims[2]=states.size();

	omega=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);


	for(p=mxGetPr(omega);it_state!=states.cend();it_state++)
	{
		memcpy(p,it_state->v.head(3).data(),sizeof(double)*dims[0]*dims[1]);
		p+=sizeof(double)*dims[0]*dims[1];
	}

	matPutVariable(result,"omega",omega);
	mxDestroyArray(omega);
}

quadrotor::DState::DState(const System & sys, const State & state, Vec4 u)
{
	Eigen::Matrix<double,M,1> dstate=f(sys,state)+h(sys,state)*u;

	v=dstate.head(6);
	a=dstate.tail(6);
}

quadopt quadrotor::optProb;
#endif
