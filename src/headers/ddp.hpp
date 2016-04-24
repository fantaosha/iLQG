/*************************************************************************
    > File Name: DDP.h
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Fri 22 Apr 2016 07:15:29 AM CDT
 ************************************************************************/
#ifndef _DDP
#define _DDP
#include <cmath>
#include <Eigen/Cholesky>
#include <algorithm>
#include <type.hpp>
#include <boxQP.hpp>


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
			// parameters to regularize 
			double lambda;
			double dlambda;
			double lambdaFactor;
			double lambdaMax;
			double lambdaMin;

			// parameter for line search
			double alpha;
			double dalphaFactor;
			double alphaMin;

			// minimal reduction ration
			double reductionRatioMin;

			// exit criterion
			double tolFun;
			double tolGrad;

			// cost function
			MatMM Q;
			MatMM Qf;
			MatNN R;

			Params(MatMM const & Q_=MatMM::Identity(),   MatNN const & R_=MatNN::Identity(),   MatMM const & Qf_=50*MatMM::Identity(),
				          double lambda_=1, double lambdaMin_=1e-5, double lambdaMax_=1e10, double dlambda_=1, double lambdaFactor_=1.6,
						  double  alpha_=1, double  alphaMin_=1e-3, double  dalphaFactor_=pow(1e-3,0.1), double reductionRatioMin_=0,
						  double  tolFun_=1e-7, double  tolGrad_ =1e-4): Q(Q_), Qf(Qf_), R(R_),
																		 lambda(lambda_), lambdaMin(lambdaMin_), lambdaMax(lambdaMax_), dlambda(dlambda_), lambdaFactor(lambdaFactor_), 
																		 alpha(alpha_), dalphaFactor(dalphaFactor_), reductionRatioMin(reductionRatioMin_), 
																		 tolFun(tolFun_), tolGrad(tolGrad_) 
			{
			}
		};

protected:
		System const sys;
		double const dt;

		Params params; // parameters to control iLQG

		int num; // number of steps

		State x0; // initial state

		std::vector<State> xs; // states
		std::vector<Ref> xrefs; // reference states

		std::vector<U> us; // control inputs 
		U umin;
		U umax;

		double J0;
		Vec2 dJ;
		
		double getGnorm(std::vector<VecN> const & kus, std::vector<VecN> const & us);
	public:               
		DDP(System const & sys, double const & dt);
		bool init(State const & xs0, std::vector<U> const & us0, std::vector<Ref> const & xrefs, Params const & params, U const & umin, U const & umax, int const & num);
		void iterate(int const & itr_max, std::vector<U> & us);
		void forwards(std::vector<MatNM> const & Ks, std::vector<VecN> const & kus, double const & alpha, std::vector<State> & xns, std::vector<U> & uns, double & Jn);
		int backwards(std::vector<MatNM> &Ks, std::vector<VecN> &kus, double const & lambda);
};

template <typename Robot> double DDP<Robot>::getGnorm(std::vector<VecN> const & kus, std::vector<VecN> const & us0)
{
	assert(kus.size()==num);
	assert(us.size()==num);

	double sum=0;
	for(int i=0;i<num;i++)
	{
		sum+=(kus[i].array().abs()/(us[i].array().abs()+1)).maxCoeff();
	}

	return sum/num;
}

template <typename Robot> DDP<Robot>::DDP(System const & sys_, double const & dt_):sys(sys_), dt(dt_), num(0)
{
}

template <typename Robot> bool DDP<Robot>::init(State const & x0_, std::vector<U> const & us_, std::vector<Ref> const & xrefs_, Params const & params_, U const & umin_, U const & umax_, int const & num_)
{
	if(xrefs_.size()<num_+1)
	{
		std::cout<<"ERROR: Not enough references."<<std::endl;
		return false;
	}

	if(us_.size()<num)
	{
		std::cout<<"ERROR: Not enough control inputs."<<std::endl;
		return false;
	}

	params=params_;
	num=num_;
	x0=x0_;

	us.reserve(num);
	xs.reserve(num+1);
	xrefs.reserve(num+1);
	us.resize(num,VecN::Zero());
	xs.resize(num+1,x0);
	xrefs.resize(num+1,x0);
	
	umin=umin_;
	umax=umax_;

	for(int i=0;i<num;i++)
	{
		us[i]=clamp<N>(us_[i],umin,umax);
		xrefs[i]=xrefs_[i];
	}
	xrefs[num]=xrefs_[num];

	xs[0]=x0;

	VecM error;

	for(int i=0;i<num;i++)
	{
		error=Robot::State::diff(xs[i], xrefs[i]);
		J0+=Robot::L(params.Q,params.R,error,us[i])*dt;

		DState dstate(sys,xs[i],us[i]);
		xs[i+1]=xs[i].update(dstate,dt);
	}
	
	error=Robot::State::diff(xs[num], xrefs[num]);
	J0+=Robot::L(params.Qf,MatNN::Zero(),error,VecN::Zero());


	std::cout<<"=========================================="<<std::endl;
	std::cout<<"J0: "<<J0<<std::endl;
	std::cout<<"=========================================="<<std::endl;
	return true;
}

template<typename Robot> void DDP<Robot>::iterate(int const & itr_max, std::vector<U> & us0)
{
	std::vector<MatNM> Ks;
	std::vector<VecN> kus;

	Ks.reserve(num);
	Ks.resize(num, MatNM::Zero());
	kus.reserve(num);
	kus.resize(num, VecN::Zero());

	double lambda=params.lambda;
	double dlambda=params.dlambda;

	for(int i=0;i<itr_max;i++)
	{
		std::cout<<"========================================================================"<<std::endl;
		std::cout<<"Iteration # "<<i<<std::endl;
		std::cout<<"------------------------------------------------------------------------"<<std::endl;
		// backward pass
		bool backPassDone=false;
		while(!backPassDone)
		{
			int result=backwards(Ks,kus,lambda);

			if(result>=0)
			{
				dlambda=std::max(dlambda*params.lambdaFactor, params.lambdaFactor);
				lambda=std::max(lambda*dlambda, params.lambdaMin);

				if(lambda>params.lambdaMax)
					break;

				continue;
			}
			
			backPassDone=true;
		}
	
		double gnorm=getGnorm(kus,us);

		if(gnorm<params.tolGrad && lambda<1e-5)
		{
			dlambda=std::min(dlambda/params.lambdaFactor, 1.0/params.lambdaFactor);
			lambda=lambda*dlambda*(lambda>params.lambdaMin);
#ifdef PRINT
			std::cout<<"SUCCESS: gradient norm = "<<gnorm" < tolGrad"<<std::endl
#endif

			break;
		}

		// forward pass
		bool fwdPassDone=false;
		std::vector<State> xns;
		std::vector<U> uns;
		double Jn;
		double actual;
		double expected;

		xns.reserve(num+1);
		uns.reserve(num);
		xns.resize(num+1,x0);
		uns.resize(num,VecN::Zero());

		if(backPassDone)
		{
			double alpha=params.alpha;

			while(alpha>params.alphaMin)
			{
				forwards(Ks, kus, alpha, xns, uns, Jn);
				actual=J0-Jn;
				expected=-alpha*dJ(0)-alpha*alpha*dJ(1);
				double reductionRatio=-1;

				if(expected>0)
					reductionRatio=actual/expected;
				else
					std::cout<<"WARNING: non-positive expected reduction: should not occur"<<std::endl;

				if(reductionRatio>params.reductionRatioMin)
					fwdPassDone=true;
				break;

				alpha*=params.dalphaFactor;
			}
		}

		std::cout<<"--------------------------------------------"<<std::endl;
		std::cout<<"Results"<<std::endl;
		std::cout<<"--------------------------------------------"<<std::endl;
		if(fwdPassDone)
		{
			dlambda=std::min(dlambda/params.lambdaFactor, 1.0/params.lambdaFactor);
			lambda=lambda*dlambda*(lambda>params.lambdaMin);
			
			std::cout<<"Improved"<<std::endl;
			std::cout<<"lambda: "<<lambda<<std::endl;
			std::cout<<"dlambda: "<<dlambda<<std::endl;
			std::cout<<"Jn: "<<Jn<<std::endl;
			std::cout<<"final state: "<<xns[num].x.transpose()<<std::endl;

			xs=xns;
			us=uns;
			J0=Jn;

			if(actual<params.tolFun)
			{
#ifdef PRINT
				std::cout<<"SUCCESS: cost change = "<<actual<<" < tolFun"<<std::endl;
#endif
				break;
			}
		}
		else
		{
			dlambda=std::max(dlambda*params.lambdaFactor, params.lambdaFactor);
			lambda=std::max(lambda*dlambda, params.lambdaMin);

			std::cout<<"No step found"<<std::endl;
			std::cout<<"lambda: "<<lambda<<std::endl;
			std::cout<<"dlambda: "<<dlambda<<std::endl;
			std::cout<<"Jn: "<<Jn<<std::endl;
			
			if (lambda>params.lambdaMax)
				break;
		}
		std::cout<<"========================================================================"<<std::endl<<std::endl;
	}

	if(us0.capacity()<num)
	{
		us0.reserve(num);
		us0.resize(num,VecN::Zero());
	}

	for(int i=0;i<num;i++)
		us0[i]=us[i];
}

template<typename Robot> void DDP<Robot>::forwards(std::vector<MatNM> const &  Ks, std::vector<VecN> const & kus,  double const & alpha,
												         std::vector<State> & xns,          std::vector<U> & uns,        double &    Jn)
{
	xns.reserve(num+1);
	uns.reserve(num);
	xns.resize(num+1,x0);
	uns.resize(num,VecN::Zero());

	Jn=0;

	xns[0]=x0;

	VecM dx, err_x;

	for(int i=0;i<num;i++)
	{
		dx=Robot::State::diff(xns[i],xs[i]);
		uns[i]=clamp<N>(us[i]+kus[i]*alpha+Ks[i]*dx,umin,umax);

		err_x=Robot::State::diff(xns[i],xrefs[i]);
		Jn+=Robot::L(params.Q, params.R, err_x, uns[i])*dt;

		DState dstate(sys,xns[i],uns[i]);
		xns[i+1]=xns[i].update(dstate,dt);
	}

	err_x=Robot::State::diff(xns[num],xrefs[num]);
	Jn+=Robot::L(params.Qf, MatNN::Zero(), err_x, VecN::Zero());
}

template<typename Robot> int DDP<Robot>::backwards(std::vector<MatNM> & Ks, std::vector<VecN> & kus, double const & lambda)
{
	// Initialization
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

	MatMM const & Q=params.Q;
	MatNN const & R=params.R;
	MatMM const & Qf=params.Qf;

	Ks.reserve(num);
	Ks.resize(num, MatNM::Zero());
	kus.reserve(num);
	kus.resize(num,U::Zero());
	
	dJ=Vec2::Zero();
	
	// Start back pass
	error=Robot::State::diff(xs[num],xrefs[num]);
	
	Vx=Robot::Lx(Qf,error);
	Vxx=Robot::Lxx(Qf,error);

	Eigen::LLT<MatNN> llt;
	
	Mat Hfree;
	int result;
	Eigen::Array<int,N,1> free;

	for(int i=num-1;i>=0;i--)
	{
		error=Robot::State::diff(xs[i],xrefs[i]);

		Lx=Robot::Lx(Q,error)*dt;
		Lxx=Robot::Lxx(Q,error)*dt;
		Lu=Robot::Lu(R,us[i])*dt;
		Luu=Robot::Luu(R,us[i])*dt;

		Robot::linearize(sys,dt,xs[i],us[i],A,B);

		At=A.transpose();
		Bt=B.transpose();

		Qx=Lx + At*Vx;
		Qu=Lu + Bt*Vx;
		Quu=Luu + Bt*Vxx*B;
		Qux=Bt*Vxx*A;
		Qxu=Qux.transpose();
		Qxx=Lxx + At*Vxx*A;

		Quum=Quu+lambda*MatNN::Identity();

		llt.compute(Quum);
		if(llt.info()!=Eigen::Success)
			return i;
		
		U dumin=umin-us[i];
		U dumax=umax-us[i];

		VecN ku;
		boxQP<N>(Quum, Qu, dumin, dumax, kus[std::min(i+1,num-1)], ku, Hfree, free, result);

		if(result<1)
			return i;
		
		MatNM K=MatNM::Zero();
		if(free.sum()>0)
		{
			int rows=free.sum();
			Mat Quuf=Hfree.transpose()*Hfree;
			Mat Quxf=Mat::Zero(rows,M);

			for(int i=0,j=0;i<N && j<rows;i++)
				if(free(i))
					Quxf.row(j++)=Qux.row(i);

			Mat Lfree=-Quuf.ldlt().solve(Quxf);

			for(int i=0, j=0;i<N && j<rows; i++)
				if(free(i))
					K.row(i)=Lfree.row(j++);
		}
		
		MatMN Kt=K.transpose();
		Eigen::Matrix<double,1,N> kut=ku.transpose();

		dJ+=(Vec2()<<kut*Qu, 0.5*kut*Quu*ku).finished();


		Vxx=Qxx+Kt*Qux+Qxu*K+Kt*Quu*K;
		Vx=Qx+Kt*Qu+Qxu*ku+Kt*Quu*ku;

		kus[i]=ku;
		Ks[i]=K;
	}

//	std::cout<<"--------------------------------------------"<<std::endl;
//	std::cout<<"Backwards"<<std::endl;
//	std::cout<<"--------------------------------------------"<<std::endl;
//	for(int i=0;i<num;i++)
//		std::cout<<"ku["<<i+1<<"]: "<<kus[i].transpose()<<std::endl;

	return -1;
}

#endif
