/*************************************************************************
    > File Name: ilqg.cpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Wed 20 Apr 2016 03:54:48 AM CDT
 ************************************************************************/
#include <Eigen/Dense>
#include <liegroup.hpp>

#include <quadrotor.hpp>
#include <iostream>

#include <ilqg.hpp>
#include <simulator.hpp>

#include <string>

//#define TRACK

int main()
{
#ifdef TRACK
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

	std::list<typename quadrotor::State> list_ref;

	for(int i=0;i<N;i+=sN)
	{
		Eigen::Matrix<double,3,3> R;
		Eigen::Matrix<double,3,1> w,xq,vq;

		memcpy(R.data(),pR,sizeof(double)*R.rows()*R.cols());
		memcpy(w.data(),pw,sizeof(double)*w.rows()*w.cols());
		memcpy(xq.data(),pxq,sizeof(double)*xq.rows()*xq.cols());
		memcpy(vq.data(),pvq,sizeof(double)*vq.rows()*vq.cols());

		list_ref.emplace_back(R,xq,w,R.transpose()*vq);

		pR+=sizeof(double)*R.rows()*R.cols()*sN;
		pw+=sizeof(double)*w.rows()*w.cols()*sN;
		pxq+=sizeof(double)*xq.rows()*xq.cols()*sN;
		pvq+=sizeof(double)*vq.rows()*vq.cols()*sN;
	}

	sN=10;

	std::list<typename quadrotor::U> list_u0;

	for(int i=0;i<N;i+=sN)
	{
		Eigen::Matrix<double,4,1> u;
		memcpy(u.data(),pU,sizeof(double)*u.rows()*u.cols());
		list_u0.push_back(u.cwiseProduct(u));
		pU+=sizeof(double)*u.rows()*u.cols()*sN;
	}
#else
	quadrotor::State state_ref;
	state_ref.g.block(0,0,3,3)=SO3::exp((Vec3()<<1,2,1).finished());
	state_ref.g.block(0,3,3,1)=(Vec3()<<0,0,0).finished();

	std::list<quadrotor::State> list_ref(4000,state_ref);
	std::list<quadrotor::U> list_u0(4000, Vec4::Zero());
#endif

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
	Mat12 Mf=5*Mat12::Identity();
//	Mf.block(0,0,3,3)*=10;
//	Mf.block(3,3,3,3)*=40;
//	Mf.block(6,6,6,6)=Mf.block(0,0,6,6)/8;
	Mat12 M=Mf/8;
	Mat4 R=Mat4::Identity()*2;
	iLQG<quadrotor>::Params params(M,R,Mf);

	// Set up initial state
	quadrotor::State state0=list_ref.front();

	state0.g.block(0,3,3,1)-=(Vec3::Random()).normalized()*5;
	state0.g.block(0,0,3,3)*=SO3::exp(Vec3::Random().normalized()*1);
	state0.v.head(3)-=Vec3::Random().normalized()*0;
	state0.v.tail(3)-=Vec3::Random().normalized()*0;
	// Set up simulator
	double dt=0.01;
	Sim<quadrotor> sim(sys,dt);
	sim.init(state0);

	Vec4 umin=-6*Vec4::Ones();
	Vec4 umax=Vec4::Ones()*6;
	
	iLQG<quadrotor> ilqg(sys,dt);

	double ts=0.02;
	double T=36;
	double Tp=1.5;
	size_t SN=size_t(ts/dt+0.5);
	size_t ND=size_t(Tp/dt+0.5)+1;

	ilqg.init(state0, list_u0, list_ref, params, umin, umax, 500);
	ilqg.evaluate(-1,200,list_u0);
/*************************************************************************
	timespec T_start, T_end;

	for(double t=1e-5; t<=T; t+=ts )
	{
		clock_gettime(CLOCK_MONOTONIC,&T_start);

		ilqg.init(state0, list_u0, list_ref, params, umin, umax, ND);
		ilqg.evaluate(-1,20,list_u0);

		clock_gettime(CLOCK_MONOTONIC,&T_end);
		std::cout<<"time consumed is "<<(T_end.tv_sec-T_start.tv_sec)+(T_end.tv_nsec-T_start.tv_nsec)/1000000000.0<<"s"<<std::endl;


		for(size_t n=0;n<SN;n++)
		{
			Vec4 ut=list_u0.front();
			sim.update(ut);

			list_u0.pop_front();
			list_ref.pop_front();
		}

		quadrotor::State state=sim.get_state();
		quadrotor::State state_ref=list_ref.front();

		std::cout<<quadrotor::State::diff(state,state_ref).transpose()<<std::endl;
	}
*************************************************************************/
}
