/*************************************************************************
    > File Name: ilqg.hpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Wed 20 Apr 2016 03:48:14 PM CDT
 ************************************************************************/
#ifndef ILQG
#define ILQG

template<typename Robot> class iLQG
{
	public:
		typedef typename Robot::System System;
		typedef typename Robot::State State;
		typedef typename Robot::DState DState;

		typedef Robot::M M;
		typedef Robot::N N;

		typedef Eigen::Matrix<double,N,1> U;
		typedef Eigen::Matrix<double,M,1> VecM;
		typedef Eigen::Matrix<double,N,1> VecN;
		typedef Eigen::Matrix<double,M,M> MatMM;
		typedef Eigen::Matrix<double,M,N> MatMN;
		typedef Eigen::Matrix<double,N,M> MatNM;
		typedef Eigen::Matrix<double,N,N> MatNN;

	protected:
		static void iLQG_back(           System const & sys,                       double      dt, 
				               std::list<State> const & list_state0, std::list<U> const & list_u0, std::list<Ref> const & list_ref, 
							              MatMM const &           M,        MatNN const &       R,          MatMM const &       Mf, 
							                  U const &        umin,            U const &    umax,         double const &     alpha, 
										        MatNM & list_K,         std::list<VecN> & list_k)
		{
			list_K.clear();
			list_k.clear();

			std::list<State>::const_reverse_iterator rit_state=list_state0.crbegin();
			std::list<Ref>::const_reverse_iterator rit_ref=list_ref.crbegin();
			std::list<U>::const_reverse_iterator rit_u=list_u0.crbegin();

			VecM dg=quadrotor::diff(*rit_state,*rit_ref);		

			MatMM Vxx=quadcost::Lxx(Mf,dg);
			VecM Vx=quadcost::Lx(Mf,dg);

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

				Vec4 Lu=Robot::Lu(R,u)*dt;
				Mat4 Luu=Robot::Luu(R,u)*dt;
				Vec12 Lx=Robot::Lx(M,dg)*dt;
				Mat12 Lxx=Robot::Lxx(M,dg)*dt;

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
}

#endif
