/*************************************************************************
    > File Name: LQR.hpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Fri 22 Apr 2016 12:20:55 PM CDT
 ************************************************************************/
#ifndef _LQR  
#define _LQR
#include <algorithm>
#include <type.hpp>


template<typename Robot> class LQR
{
	public:
		typedef typename Robot::System System;
		typedef typename Robot::State State;
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
			MatMM Q;
			MatMM Qf;
			MatNN R;

			Params(MatMM const & Q_=MatMM::Identity(), MatNN const & R_=MatNN::Identity(), MatMM const & Qf_=50*MatMM::Identity()): Q(Q_), Qf(Qf_), R(R_)
			{
			}
		};

	public:               
		bool static solve (System const & sys, double const & dt, size_t const & num, Params const & params, std::vector<U> const us, std::vector<Ref> const refs, std::vector<MatNM> & Ks)
		{
			if(us.size()<num || refs.size()<(num+1))
			{
				std::cout<<"Not enough control inputs and states"<<std::endl;
				return false;
			}

			Ks.clear();
			Ks.reserve(num);

			Ref ref;
			U u;

			MatMM A;
			MatMN B;

			MatMM const & Q=params.Q;
			MatMM const & Qf=params.Qf;
			MatNN const & R=params.R;

			MatMM P=Qf;

			for(int i=num-1;i>=0;i--)
			{
				ref=refs[i];
				u=us[i];

				Robot::linearize(sys, dt, ref, u, A, B);
				MatMM At=A.transpose();
				MatNM Bt=B.transpose();

				Eigen::Matrix<double,4,12> K=(R+Bt*P*B).inverse()*Bt*P*A;
				P=K.transpose()*R*K+Q+(A-B*K).transpose()*P*(A-B*K);

				Ks[i]=K;
			}

			return true;
		}
};

#endif
