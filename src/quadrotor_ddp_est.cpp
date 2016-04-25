/*************************************************************************
    > File Name: quadrotor_ddp.cpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Fri 22 Apr 2016 09:27:35 AM CDT
 ************************************************************************/
#include <Eigen/Dense>
#include <liegroup.hpp>

#include <quadrotor.hpp>
#include <iostream>

#include <ddp.hpp>
#include <simulator.hpp>

#include <string>
#define TRACK

int main()
{
	char file[]="/media/fantaosha/Documents/JHU/Summer 2015/quad_rotor_traj/traj_full_12s_1_small.mat";

	mxArray *mxR;
	mxArray *mxw;
	mxArray *mxxq;
	mxArray *mxvq;
	mxArray *mxU;

	mxR=matGetVariable(matOpen(file,"r"),"state_R");
	mxw=matGetVariable(matOpen(file,"r"),"state_w");
	mxxq=matGetVariable(matOpen(file,"r"),"state_xq");
	mxvq=matGetVariable(matOpen(file,"r"),"state_vq");

	mxU=matGetVariable(matOpen(file,"r"),"U");

	size_t N=mxGetN(mxU);
	size_t sN=10;

	void *pR=mxGetPr(mxR);
	void *pw=mxGetPr(mxw);
	void *pxq=mxGetPr(mxxq);
	void *pvq=mxGetPr(mxvq);
	void *pU=mxGetPr(mxU);

	std::vector<typename quadrotor::State> xrefs0;
	xrefs0.reserve(N/sN+1);

	for(int i=0;i<N;i+=sN)
	{
		Eigen::Matrix<double,3,3> R;
		Eigen::Matrix<double,3,1> w,xq,vq;

		memcpy(R.data(),pR,sizeof(double)*R.rows()*R.cols());
//		R=(Vec3()<<1,-1,-1).finished().asDiagonal();
//		w<<0,0,0;
		memcpy(w.data(),pw,sizeof(double)*w.rows()*w.cols());
		memcpy(xq.data(),pxq,sizeof(double)*xq.rows()*xq.cols());
		memcpy(vq.data(),pvq,sizeof(double)*vq.rows()*vq.cols());

		xrefs0.emplace_back(R,xq,w,R.transpose()*vq);

		pR+=sizeof(double)*R.rows()*R.cols()*sN;
		pw+=sizeof(double)*w.rows()*w.cols()*sN;
		pxq+=sizeof(double)*xq.rows()*xq.cols()*sN;
		pvq+=sizeof(double)*vq.rows()*vq.cols()*sN;
	}

	sN=10;

	std::vector<typename quadrotor::U> us0;
	us0.reserve(N/sN+1);

	for(int i=0;i<N;i+=sN)
	{
		Eigen::Matrix<double,4,1> u;
		memcpy(u.data(),pU,sizeof(double)*u.rows()*u.cols());
		us0.push_back(u.cwiseProduct(u));
		pU+=sizeof(double)*u.rows()*u.cols()*sN;
	}


	std::srand((unsigned int) time(0));
	// Set up quadrotor parameters
	double m=0.6;
	double Ix=8*1e-3; // moment of inertia in X- Y- and Z-axis
	double Iy=7.5*1e-3;
	double Iz=13.5*1e-3;
	double d=0.2; // displacement of rotors
	double kt=0.6;
	double km=0.15;

	quadrotor::System sys(Ix,Iy,Iz,m,d,km,kt);
	
	// Set up cost function
	Mat12 Mf=5*Mat12::Identity()/10;
	Mf.block(0,0,3,3)*=1;
	Mf.block(3,3,3,3)*=6;
//	Mf.block(6,6,6,6)=Mf.block(0,0,6,6);
	Mat12 M=Mf/2;
	Mat4 R=Mat4::Identity()/50;
	DDP<quadrotor>::Params params(M,R,Mf);
	double dt=0.01;

	// Set up initial state

	int NUM=10;
	int succeed=0;
//#pragma omp parallel for reduction(+:succeed) schedule(dynamic)  
	for(int n=0;n<NUM;n++)
	{
			std::vector<typename quadrotor::State> xrefs=xrefs0;
			std::vector<typename quadrotor::U> us=us0;
			quadrotor::State x0=xrefs[0];

			x0.g.block(0,3,3,1)-=(Vec3::Random()).normalized()*30;
			x0.g.block(0,0,3,3)*=SO3::exp(Vec3::Random().normalized()*3);
		//	x0.g.block(0,0,3,3)*=SO3::exp((Vec3()<<3.14,0,0).finished());
			x0.v.head(3)-=Vec3::Random().normalized()*0;
			x0.v.tail(3)-=Vec3::Random().normalized()*0;
			// Set up simulator
			
			Sim<quadrotor> sim(sys,dt);
			sim.init(x0,4000);

			Vec4 umin=-Vec4::Ones()*0;
			Vec4 umax=Vec4::Ones()*6;
			
			DDP<quadrotor> ddp(sys,dt);

			double ts=0.02;
			double T=36;
			double Tp=1.5;
			size_t SN=size_t(ts/dt+0.5);
			size_t ND=size_t(Tp/dt+0.5)+1;

#pragma omp critical
			{
//				std::cout<<quadrotor::State::diff(sim.get_state(),xrefs[0]).transpose()<<std::endl;
//				std::cout<<xrefs.size()<<std::endl;
//				std::cout<<us.size()<<std::endl;
			}

			int sn=2;
			int itr_max=15;
			for(int i=0;i<4000;i+=sn)
			{
				timespec T_start, T_end;
				clock_gettime(CLOCK_MONOTONIC,&T_start);
				ddp.init(sim.get_state(), us, xrefs, params, umin, umax, 150);
				ddp.iterate(itr_max,us);
				clock_gettime(CLOCK_MONOTONIC,&T_end);
//				std::cout<<"time consumed is "<<(T_end.tv_sec-T_start.tv_sec)+(T_end.tv_nsec-T_start.tv_nsec)/1000000000.0<<"s"<<std::endl;
			
				for(int j=0;j<sn;j++)
				{
					sim.update(us.front());
					xrefs.erase(xrefs.begin());
					us.erase(us.begin());
				}

				Vec12 error=quadrotor::State::diff(sim.get_state(),xrefs[0]);

				if(error.norm()<5)
				{
					succeed++;
					break;
				}
			}
#pragma omp critical
			{
				std::cout<<quadrotor::State::diff(sim.get_state(),xrefs[0]).transpose()<<std::endl;
			}
	}
	
	std::cout<<succeed<<" "<<NUM<<std::endl;
}
