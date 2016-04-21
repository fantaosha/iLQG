/*************************************************************************
    > File Name: ilqg.hpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Wed 20 Apr 2016 03:48:14 PM CDT
 ************************************************************************/
#ifndef _ILQG
#define _ILQG
#include <liegroup.hpp>
#include <snOPT.hpp>
#include <Eigen/Cholesky>
#include <algorithm>

template<typename Robot> class iLQG
{
	public:
		typedef typename Robot::System System;
		typedef typename Robot::State State;
		typedef typename Robot::DState DState;
		typedef typename Robot::State Ref;

		static const size_t M=Robot::M;
		static const size_t N=Robot::N;
		
		typedef Eigen::Matrix<double,M,1> VecM;
		typedef Eigen::Matrix<double,N,1> VecN;
		typedef Eigen::Matrix<double,N,1> U;
		typedef Eigen::Matrix<double,M,M> MatMM;
		typedef Eigen::Matrix<double,M,N> MatMN;
		typedef Eigen::Matrix<double,N,M> MatNM;
		typedef Eigen::Matrix<double,N,N> MatNN;

		struct Params
		{
			double mu;

			double dmu0;

			double mumax;
			double mumin;

			MatMM M;
			MatMM Mf;
			MatNN R;

			Params(MatMM const & M_=MatMM::Identity(), MatNN const & R_=MatNN::Identity(), MatMM const & Mf_=50*MatMM::Identity(),
				   double mu_=1,  double dmu0_=1.6, double mumin_=1e-5, double mumax_=1e5): mu(mu_), mumin(mumin_), mumax(mumax_), dmu0(dmu0_), M(M_), Mf(Mf_), R(R_)
			{
			}
		};

	protected:
		typedef	snOPT<N> snopt; 

protected:
		System const sys;
		double const dt;

		Params params;

		size_t Num;

		State state0;
		
		VecN umax,umin;

		std::list<U> list_u0;
		std::list<U> list_u1;
		std::list<Ref> list_ref;
		std::list<State> list_state0;

		double J0;

		int box_opt(MatNN const & Quu, VecN const & Qu, VecN const & dumin, VecN const & dumax, VecN & ku, Eigen::Matrix<int,N,1> & ku_state);
		
		bool box(VecN &x, Eigen::Matrix<int,N,1> & xstate, VecN const & xmin, VecN const & xmax);

	public:
		iLQG(System const & sys, double const & dt);
		bool init(State const & state0, std::list<U> const & list_u0, std::list<Ref> const & list_ref, Params const & params, VecN const & umin, VecN const & umax, size_t const & Num);
		void evaluate(double const & tolerance, size_t const & itr_max, std::list<U> & list_u);
		double forwards(std::list<MatNM> const & list_K, std::list<VecN> const & list_ku);
		void backwards(std::list<MatNM> &list_K, std::list<VecN> &list_ku);
};

template <typename Robot> iLQG<Robot>::iLQG(System const & sys_, double const & dt_):sys(sys_), dt(dt_), Num(0)
{
	snopt::init();
}

template <typename Robot> bool iLQG<Robot>::init(State const & state0_, std::list<U> const & list_u0_, std::list<Ref> const & list_ref_, Params const & params_, VecN const & umin_, VecN const & umax_, size_t const & Num_)
{
	if(list_ref_.size()<Num_)
	{
		std::cout<<"ERROR: Not enough references."<<std::endl;
		return false;
	}

	if(list_u0_.size()<Num_-1)
	{
		std::cout<<"ERROR: Not enough control inputs."<<std::endl;
		return false;
	}

	params=params_;
	Num=Num_;
	state0=state0_;
	umin=umin_;
	umax=umax_;

	list_u0.clear();
	list_u1.clear();
	list_ref.clear();
	list_state0.clear();

	double J0=0;

	typename std::list<U>::const_iterator itr_u=list_u0_.begin();
	typename std::list<Ref>::const_iterator itr_ref=list_ref_.begin();

	State state=state0;
	list_state0.push_back(state);

	VecM error;

	for(int i=1;i<Num;i++)
	{
		error=Robot::State::diff(state, *itr_ref);
		J0+=Robot::L(params.M,params.R,error,*itr_u)*dt;

		DState dstate(sys,state,*itr_u);
		state=state.update(dstate,dt);

		list_state0.push_back(state);
		list_u0.push_back(*itr_u);
		list_ref.push_back(*itr_ref);

		itr_u++;
		itr_ref++;
	}
	
	list_ref.push_back(*itr_ref);
	error=Robot::State::diff(state, *itr_ref);
	J0+=Robot::L(params.Mf,MatNN::Zero(),error,VecN::Zero());


	std::cout<<"=========================================="<<std::endl;
	std::cout<<J0<<std::endl;
	std::cout<<state.x<<std::endl;
	std::cout<<"=========================================="<<std::endl;
	return true;
}

template<typename Robot> void iLQG<Robot>::evaluate(double const & tolerance, size_t const & itr_max, std::list<U> & list_u)
{
	std::list<MatNM> list_K;
	std::list<VecN> list_ku;

	for(int i=0; i<itr_max; i++)
	{
		backwards(list_K, list_ku);
		double J1=forwards(list_K,list_ku);

		if(fabs(J1-J0)/J0<tolerance)
			break;

		J0=J1;
	}

	typename std::list<U>::const_iterator itr_u0=list_u0.begin();
	typename std::list<U>::iterator itr_u=list_u.begin();

	while(itr_u0!=list_u0.end())
	{
		*itr_u=*itr_u0;

		itr_u0++;
		itr_u++;
	}
}

template<typename Robot> double iLQG<Robot>::forwards(std::list<MatNM> const & list_K, std::list<VecN> const & list_ku)
{
	std::list<State> list_state;

	typename std::list<U>::const_iterator itr_u0=list_u0.begin();
	typename std::list<State>::const_iterator itr_state0=list_state0.begin();
	typename std::list<State>::const_iterator itr_ref=list_ref.begin();
	typename std::list<U>::const_iterator itr_ku=list_ku.begin();
	typename std::list<MatNM >::const_iterator itr_K=list_K.begin();
	
	State state=list_state0.front();

	State state0;
	Ref ref;
	VecN u0;
	VecN u1;

	VecN ku;
	MatNM K;

	VecM error;
	VecM dg;

	double J1=0;	
	std::cout<<"============================================"<<std::endl;
	for(int i=1;i<Num;i++)
	{
		list_state.push_back(state);

		state0=*itr_state0;
		ref=*itr_ref;
		u0=*itr_u0;
	

		ku=*itr_ku;
		K=*itr_K;


		dg=Robot::State::diff(state,state0);
		u1=u0+ku+K*dg;
		
		std::cout<<u1.transpose()<<std::endl;
/*************************************************************************
		for(int i=0;i<4;i++)
			if(u1(i)>umax(i))
				u1(i)=umax(i);
			else
				if(u1(i)<umin(i))
					u1(i)=umin(i);
*************************************************************************/

		error=Robot::State::diff(state, ref);
		J1+=Robot::L(params.M,params.R,error,u1)*dt;
//		std::cout<<J1<<" ";
//		if(J1>1000)
//		{
//			std::cout<<"***********************************"<<std::endl;
//			std::cout<<state.v.transpose()<<std::endl;
//			std::cout<<ref.v.transpose()<<std::endl;
//			std::cout<<"***********************************"<<std::endl;
//		}

		DState dstate(sys,state,u1);
		state=state.update(dstate,dt);
		
		list_u1.push_back(u1);

		itr_u0++;
		itr_state0++;
		itr_ref++;
		itr_ku++;
		itr_K++;

	}

	error=Robot::State::diff(state, ref);
	J1+=Robot::L(params.Mf,MatNN::Zero(),error,U::Zero());

	list_state.push_back(state);

	list_state0=list_state;
	list_u0=list_u1;

	std::cout<<"J1: "<<J1<<std::endl;
	std::cout<<state.x<<std::endl;
	std::cout<<"============================================"<<std::endl;
//	std::cout<<state.g<<std::endl;

	return J1;
}

template<typename Robot> void iLQG<Robot>::backwards(std::list<MatNM> &list_K, std::list<VecN> &list_ku)
{
	VecM Qx;
	VecN Qu;
	MatMM Qxx;
	MatMN Qxu;
	MatNM Qux;
	MatNN Quu;
	MatNN Quum;

	VecM Lx;
	VecN Lu;
	MatMM Lxx;
	MatNN Luu;

	MatMM Vxx;
	VecM Vx;

	MatMM A;
	MatMM At;
	MatMN B;
	MatNM Bt;

	VecM error;

	MatMM const & M=params.M;
	MatNN const & R=params.R;
	MatMM const & Mf=params.Mf;
	
	typename std::list<State>::const_reverse_iterator rit_state=list_state0.crbegin();
	typename std::list<Ref>::const_reverse_iterator rit_ref=list_ref.crbegin();
	typename std::list<U>::const_reverse_iterator rit_u=list_u0.crbegin();

	error=Robot::State::diff(*rit_state,*rit_ref);
	Vx=Robot::Lx(Mf,error);
	Vxx=Robot::Lxx(Mf,error);

	rit_state++;
	rit_ref++;

	list_K.clear();
	list_ku.clear();

	Eigen::LLT<MatNN> llt;
	while(rit_state!=list_state0.crend())
	{
		State const &state=*rit_state;
		U const & u=*rit_u;
		Ref const & ref=*rit_ref;

		error=Robot::State::diff(*rit_state,*rit_ref);
		Lx=Robot::Lx(M,error)*dt;
		Lu=Robot::Lu(R,u)*dt;
		Lxx=Robot::Lxx(M,error)*dt;
		Luu=Robot::Luu(R,u)*dt;

		Robot::linearize(sys, dt, state, u, A, B);
	
//		std::cout<<A<<std::endl;
//		std::cout<<B<<std::endl;
		At=A.transpose();
		Bt=B.transpose();


		Qx=Lx + At*Vx;
		Qu=Lu + Bt*Vx;

		Qxx=Lxx + At*Vxx*A;
		Quu=Luu + Bt*Vxx*B;
		Qux=Bt*Vxx*A;
		Qxu=Qux.transpose();

//		std::cout<<Qx.transpose()<<std::endl;
//		std::cout<<Qu.transpose()<<std::endl;
//		std::cout<<Qxx<<std::endl;
//		std::cout<<Qux<<std::endl;
//		std::cout<<Quu<<std::endl;

		double mu=params.mu;
		double dmu=1;


		while(1)
		{
			Quum=Quu+mu*MatNN::Identity();
			llt.compute(Quum);
			
			if(llt.info()==Eigen::Success)
			{
				dmu=std::min(1/params.dmu0,dmu/params.dmu0);
				mu=std::max(mu*dmu,params.mumin);

				break;
			}

			dmu=std::max(params.dmu0, dmu*params.dmu0);
			mu=std::max(params.mumin,mu*dmu);
		}
		
		if (mu>params.mumax)
			break;
		
//		std::cout<<"mu: "<<mu<<std::endl;

		VecN ku=-llt.solve(Qu);
		MatNM K=-llt.solve(Qux); 


/*************************************************************************
		VecN dumax=umax-u;
		VecN dumin=umin-u;

		Eigen::Matrix<int,N,1> ku_state=Eigen::Matrix<int,N,1>::Zero();

//		std::cout<<ku.transpose()<<std::endl;
		
		if(box(ku,ku_state,dumin,dumax))
		{
			VecN ku0=ku;
			Eigen::Matrix<int,N,1> ku_state0;
			int result=box_opt(Quum, Qu, dumin, dumax, ku0, ku_state0);
			if(result>=2)
			{
				double F=ku.transpose()*(Quum*ku*0.5+Qu);
				double F0=ku0.transpose()*(Quum*ku0*0.5+Qu);

				if(F0<F)
				{
					ku_state=ku_state0;
					ku=ku0;
				}
			}
			else
			{
				ku_state=ku_state0;
				ku=ku0;
			}

			for(int i=0;i<N;i++)
				if(ku_state(i))
					K.row(i)=Eigen::Matrix<double,1,iLQG<Robot>::M>::Zero();
		}
*************************************************************************/

		MatMN Kt=K.transpose();


		Vxx=Qxx+Kt*Quu*K+Kt*Qux+Qxu*K;
		Vx=Qx+Kt*(Quu*ku+Qu)+Qxu*ku;

		std::cout<<K<<std::endl;
		std::cout<<ku.transpose()<<std::endl<<std::endl;
//		std::cout<<Vxx<<std::endl;
//		std::cout<<Vx.transpose()<<std::endl;

		list_ku.push_front(ku);
		list_K.push_front(K);

		rit_state++;
		rit_ref++;
		rit_u++;
	}
	
//	std::cout<<"size: "<<list_ku.size()<<std::endl;
}

template<typename Robot> int iLQG<Robot>::box_opt(MatNN const & Quu, VecN const & Qu, VecN const & dumin, VecN const & dumax, VecN & ku, Eigen::Matrix<int,N,1>  & ku_state)
{
	typename snopt::Params options;
	options.xlow=dumin;
	options.xupp=dumax;

	VecN x=ku;
	VecN xmul=VecN::Zero();
	Eigen::Matrix<int,N,1> xstate=Eigen::Matrix<int,N,1>::Zero();

	double F=0;
	double Fmul=0;
	int Fstate=0;

	int result=snopt::fmin(Quu,   Qu, options,
			                 x, xmul,  xstate,
							 F, Fmul,  Fstate);

	box(x,xstate,options.xlow,options.xupp);

	ku=x;
	ku_state=xstate;

	return result;
}

template<typename Robot> bool iLQG<Robot>::box(VecN &x, Eigen::Matrix<int,N,1> & xstate, VecN const & xmin, VecN const & xmax)
{
	bool bounded=false;
	xstate=Eigen::Matrix<int,N,1>::Zero();

	for(int i=0;i<Robot::N;i++)
		if(x(i)>xmax(i))
		{
			x(i)=xmax(i);
			xstate(i)=1;

			bounded=true;
		}
		else
			if(x(i)<xmin(i))
			{
				x(i)=xmin(i);
				xstate(i)=1;

				bounded=true;
			}

	return bounded;
}
#endif
