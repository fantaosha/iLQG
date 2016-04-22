/*************************************************************************
    > File Name: DDP.h
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Fri 22 Apr 2016 07:15:29 AM CDT
 ************************************************************************/
#ifndef _DDP
#define _DDP
#include <liegroup.hpp>
#include <snOPT.hpp>
#include <Eigen/Cholesky>
#include <algorithm>
#include <type.hpp>


template<typename Robot> class DDP
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
		System const sys;
		double const dt;

		Params params;

		size_t Num;

		State state0;

		std::list<U> list_u0;
		std::list<U> list_du;
		std::list<Ref> list_ref;
		std::list<State> list_state0;

		double J0;

		Vec2 dV;

		constexpr static  double dJmin=1e-10;
		const double a; // Armijo search stepsize
		constexpr static  double amin=1e-10; 

		constexpr static  double sigma=0.1;
		constexpr static  double beta=0.25; 

		constexpr static double s1=0.1;   ///< Armijo-Goldstein step-size control factor s1
		constexpr static double s2=0.5;   ///< Armijo-Goldstein step-size control factor s2
		constexpr static double b1=0.25;   ///< Armijo-Goldstein step-size control factor b1
		constexpr static double b2=2;   ///< Armijo-Goldstein step-size control factor b2

	public:               
		DDP(System const & sys, double const & dt);
		bool init(State const & state0, std::list<U> const & list_u0, std::list<Ref> const & list_ref, Params const & params, size_t const & Num);
		void iterate(double const & tolerance, size_t const & itr_max, std::list<U> & list_u);
		void forwards(std::list<MatNM> const & list_K, std::list<VecN> const & list_ku);
		void backwards(std::list<MatNM> &list_K, std::list<VecN> &list_ku);
};

template <typename Robot> DDP<Robot>::DDP(System const & sys_, double const & dt_):sys(sys_), dt(dt_), Num(0), a(2)
{
}

template <typename Robot> bool DDP<Robot>::init(State const & state0_, std::list<U> const & list_u0_, std::list<Ref> const & list_ref_, Params const & params_, size_t const & Num_)
{
	if(list_ref_.size()<Num_+1)
	{
		std::cout<<"ERROR: Not enough references."<<std::endl;
		return false;
	}

	if(list_u0_.size()<Num)
	{
		std::cout<<"ERROR: Not enough control inputs."<<std::endl;
		return false;
	}

	params=params_;
	Num=Num_;
	state0=state0_;

	list_u0.clear();
	list_du.clear();
	list_ref.clear();
	list_state0.clear();

//	alpha=1;
	dV=Vec2::Zero();

	J0=0;

	typename std::list<U>::const_iterator itr_u=list_u0_.begin();
	typename std::list<Ref>::const_iterator itr_ref=list_ref_.begin();

	State state=state0;
	list_state0.push_back(state);

	VecM error;

	for(int i=0;i<Num;i++)
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
	std::cout<<"J0: "<<J0<<std::endl;
	std::cout<<"=========================================="<<std::endl;
	return true;
}

template<typename Robot> void DDP<Robot>::iterate(double const & tolerance, size_t const & itr_max, std::list<U> & list_u)
{
	std::list<MatNM> list_K;
	std::list<VecN> list_ku;

	for(int i=0; i<itr_max; i++)
	{
		backwards(list_K, list_ku);
		forwards(list_K,list_ku);
		
		typename std::list<U>::iterator itr_u0=list_u0.begin();
		typename std::list<U>::const_iterator itr_du=list_du.cbegin();

		while(itr_du!=list_du.cend())
			*itr_u0+=*itr_du;

		J0=0;
		typename std::list<Ref>::const_iterator itr_ref=list_ref.cbegin();

		State state=state0;
		list_state0.clear();
		list_state0.push_back(state);

		VecM error;

		itr_u0=list_u0.begin();
		for(int i=0; i<Num;i++)
		{
			U u0=*(itr_u0);
			Ref ref=*itr_ref;

			error=Robot::State::diff(state, ref);
			J0+=Robot::L(params.M,params.R,error,u0)*dt;

			DState dstate(sys,state,u0);
			state=state.update(dstate,dt);

			list_state0.push_back(state);

			itr_u0++;
			itr_ref++;
		}
		
		list_ref.push_back(*itr_ref);
		error=Robot::State::diff(state, *itr_ref);
		J0+=Robot::L(params.Mf,MatNN::Zero(),error,VecN::Zero());

		std::cout<<"===================================================="<<std::endl;
		std::cout<<"Iteration "<<i<<std::endl;
		std::cout<<"J0: "<<J0<<std::endl;
		std::cout<<"===================================================="<<std::endl;
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

template<typename Robot> void DDP<Robot>::forwards(std::list<MatNM> const & list_K, std::list<VecN> const & list_ku)
{
	double a=this->a;
	double dJ=1;

	for(int count=20;count>0;count--)
	{
		typename std::list<U>::const_iterator itr_u0=list_u0.begin();
		typename std::list<State>::const_iterator itr_state0=list_state0.begin();
		typename std::list<State>::const_iterator itr_ref=list_ref.begin();
		typename std::list<U>::const_iterator itr_ku=list_ku.begin();
		typename std::list<MatNM >::const_iterator itr_K=list_K.begin();


		Ref ref;
		VecN u0;
		VecN du;
		VecN un;

		VecN ku;
		MatNM K;

		VecM error;
		VecM dg;

		double Jn=0;

		list_du.clear();
	
//		std::cout<<"****************************"<<std::endl;
//		std::cout<<"Minior Iteration "<<11-count<<std::endl;
//		std::cout<<"****************************"<<std::endl;
		State state=list_state0.front();

		for(int i=0;i<Num;i++)
		{
			state0=*itr_state0;
			ref=*itr_ref;
			u0=*itr_u0;

			ku=*itr_ku;
			K=*itr_K;

			dg=Robot::State::diff(state,state0);
			
			if(std::isnan(dg(0)))
			{
				std::cout<<"Iteration #"<<i<<std::endl;
				assert(std::isnan(dg(0)));
			}
			du=a*ku+K*dg;
			un=u0+du;

			error=Robot::State::diff(state, ref);
			Jn+=Robot::L(params.M,params.R,error,un)*dt;

			try
			{
				DState dstate(sys,state,un);
				state=state.update(dstate,dt);

				if(std::isnan(state.g(0)))
					throw 1;
			}
			catch(int e)
			{
//				std::cerr << "exception caught: " << e.what() <<std::endl; 
				std::cout << Jn<<std::endl;
				std::cout << "NaN is observed"<<std::endl;
//				throw std::runtime_error(std::string("Nan observed")); 

				return;
			}

			list_du.push_back(du);

			itr_u0++;
			itr_state0++;
			itr_ref++;
			itr_ku++;
			itr_K++;
		}

		dg=Robot::State::diff(state,state0);

		error=Robot::State::diff(state, ref);
		Jn+=Robot::L(params.Mf,MatNN::Zero(),error,U::Zero());

		double dJ=Jn-J0;
	
		if(a<amin || fabs(dJ) < dJmin)
			break;

		if(dJ>sigma*a*dV(0))
		{
			a*=beta;
			continue;
		}
		else
			break;

/*************************************************************************
		if(dJ>0)
		{
			a*=b1;
		
			if(a<1e-12)
				break;
		}
		else
		{
			double dJ0=a*dV(0)+0.5*a*a*dV(1);
			double dJm=dJ/dJ0;

			if(dJm<s1)
				a*=b1;
			else
				if (dJm>s2)
					alpha*=dalpha2;
				else
					break;
			if(a<1e-12)
				break;
		}
*************************************************************************/
	}

//	std::cout<<"============================================"<<std::endl;
//	std::cout<<"J1: "<<J1<<std::endl;
//	std::cout<<list_state.back().x.transpose()<<std::endl;
//	std::cout<<"============================================"<<std::endl;
//	std::cout<<state.g<<std::endl;
}

template<typename Robot> void DDP<Robot>::backwards(std::list<MatNM> &list_K, std::list<VecN> &list_ku)
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

	dV=Vec2::Zero();

	static const MatNN Id=MatNN::Identity();

	typename std::list<State>::const_reverse_iterator rit_state=list_state0.crbegin();
	typename std::list<Ref>::const_reverse_iterator rit_ref=list_ref.crbegin();
	typename std::list<U>::const_reverse_iterator rit_u=list_u0.crbegin();

	error=Robot::State::diff(*rit_state,*rit_ref);
	Vx=Robot::Lx(Mf,error);
	Vxx=Robot::Lxx(Mf,error);
	
	list_K.clear();
	list_ku.clear();

	Eigen::LLT<MatNN> llt;

	rit_state++;
	rit_ref++;

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

			if (mu>params.mumax)
				break;
		}
		
		if (mu>params.mumax)
			break;
		
		VecN ku=-llt.solve(Qu);
		MatNM K=-llt.solve(Qux); 

		assert(!std::isnan(ku[0]));
		assert(!std::isnan(K(0,0)));

		MatMN Kt=K.transpose();


		Vxx=Qxx+Kt*Quum*K+Kt*Qux+Qxu*K;
		Vx=Qx+Kt*(Quum*ku+Qu)+Qxu*ku;
		
		dV(0)+=ku.dot(Qu);
		dV(1)+=0.5*ku.dot(Quu*ku);

		list_ku.push_front(ku);
		list_K.push_front(K);

		rit_state++;
		rit_ref++;
		rit_u++;
	}
}

#endif
