/*************************************************************************
    > File Name: quadrotor.hpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Tue 21 Jul 2015 11:10:23 PM EDT
 ************************************************************************/
#ifndef QUADROTOR
#define QUADROTOR

#include<iostream>
#include<cmath>
#include<string>
#include<functional>

#include<Eigen/Dense>
#include<mat.h>
#include<liegroup.hpp>

struct quadrotor
{
	public:
		struct System;
		struct State;
		struct DState;

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
			static void  save(const std::vector<State> &, std::string);
			static Vec12 diff(State const & state, State const & ref);
		};
		
		struct DState // time derivative of quadrotor state
		{
			Vec6 v; // body twist
			Vec6 a; // body twist accelebration;

			DState(const System & sys, const State & state, Vec4 u); // compute body velocity
			DState(DState k1, DState k2, DState k3, DState k4); // compute average body velocity for RK4
		};

	public:
		typedef State Ref; // references used to compute cost function

		static constexpr double g=9.8; // accelebration of gravity

		static constexpr size_t M=12; // dimension of the configuration space
		static constexpr size_t N=4;  // dimension of control inputs

		typedef Eigen::Matrix<double, M,1> V;
		typedef Eigen::Matrix<double, N,1> U;


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

		static void linearize (System const & sys, double const & dt,  State const & state, U const & u, Eigen::Matrix<double,12,12> & A, Eigen::Matrix<double,12,4> & B)
		{
			Mat6 Ad=SE3::Ad(SE3::exp(-state.v*dt));
			Mat6 dexp=SE3::dexp(state.v*dt);
			
			Mat12 dgF=quadrotor::Dgf(sys,state,u);
			dgF.block(0,0,6,12)=dexp*dgF.block(0,0,6,12);

			A=Mat12::Identity()+dgF*dt;
			A.block(0,0,6,12)=Ad*A.block(0,0,6,12);
			
			Eigen::Matrix<double,M,N> duF=quadrotor::Duf(sys,state, u);
			B=duF*dt;
		}

		static double L(const Mat12 &M, const Mat4 &R, const Vec12 &dg, const Vec4 &du)
		{

			return (dg.transpose()*M*dg+du.transpose()*R*du)(0)*0.5;
		}

		static Vec12 Lx(const Mat12 &M, const Vec12 &dg)
		{
			Vec12 Lx=M*dg;
			Lx.head(6)=SE3::dexpinv(-dg.head(6)).transpose()*Lx.head(6);

			return Lx;
		}

		static Mat12 Lxx(const Mat12 &M, const Vec12 &dg)
		{
			Vec6 dg_x=dg.head(6);
			Vec6 dg_v=dg.tail(6);

			Mat6 dexpinv=SE3::dexpinv(-dg_x);
			Mat6 dexpinvT=dexpinv.transpose();

			Eigen::Matrix<double,6,36> ddexpinv=SE3::ddexpinv(-dg_x);

			Mat12 Lxx=M;
			Lxx.block(0,0,6,6)=dexpinvT*M.block(0,0,6,6)*dexpinv;
			Lxx.block(0,6,6,6)=dexpinvT*M.block(0,6,6,6);
			Lxx.block(6,0,6,6)=M.block(0,6,6,6).transpose();

/*************************************************************************
			Eigen::Matrix<double,1,6> r1=dg.transpose()*M.block(0,0,12,6);

			Mat6 DM1=-dexpinvT*(Mat6()<<r1*ddexpinv.block(0,0,6,6),
										r1*ddexpinv.block(0,6,6,6),
										r1*ddexpinv.block(0,12,6,6),
										r1*ddexpinv.block(0,18,6,6),
										r1*ddexpinv.block(0,24,6,6),
										r1*ddexpinv.block(0,30,6,6)).finished();
			
			Eigen::Matrix<double,1,6> r2=-0.5*r1*dexpinv;

			Mat6 DM2=(Mat6()<<r2*SE3::ad(SE3::e[0]),
							  r2*SE3::ad(SE3::e[1]),
							  r2*SE3::ad(SE3::e[2]),
							  r2*SE3::ad(SE3::e[3]),
							  r2*SE3::ad(SE3::e[4]),
							  r2*SE3::ad(SE3::e[5])).finished();
			
			Lxx.block(0,0,6,6)-=DM1+DM2;	
*************************************************************************/

			return Lxx;
		}

		static Vec4 Lu(const Mat4 &R, const Vec4 &du)
		{
			return R*du;
		}

		static Mat4 Luu(const Mat4 &R, const Vec4 &du)
		{
			return R;
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

void quadrotor::State::save(const std::vector<State> & states, std::string path)
{
	MATFile *result;

	mxArray *R;
	mxArray *xq;
	mxArray *omega;
	mxArray *vq;

	void *p;
	
	mwSize *dims=new mwSize[3];
	
	result=matOpen(path.c_str(),"w");

	std::vector<State>::const_iterator it_state;

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

Vec12 quadrotor::State::diff(State const & state, State const & ref)
{
	return (Vec12()<<SE3::log(SE3::inverse(ref.g)*state.g), state.v-ref.v).finished();
}

quadrotor::DState::DState(const System & sys, const State & state, Vec4 u)
{
	Eigen::Matrix<double,M,1> dstate=f(sys,state)+h(sys,state)*u;

	v=dstate.head(6);
	a=dstate.tail(6);
}

quadrotor::DState::DState(DState k1, DState k2, DState k3, DState k4)
{
	v=(k1.v+2*k2.v+2*k3.v+k4.v)/6;
	a=(k1.a+2*k2.a+2*k3.a+k4.a)/6;
}
#endif
