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

#include <string>
int main()
{
	double m=0.6;
	double Ix=8*1e-3; // moment of inertia in X- Y- and Z-axis
	double Iy=7.5*1e-3;
	double Iz=13.5*1e-3;
	double d=0.2; // displacement of rotors
	double kt=0.6;
	double km=0.15;

	quadrotor::System sys(Ix,Iy,Iz,m,d,km,kt);

	double dt=0.01;
	double T=0.2;

	quadrotor::State state0(Mat4::Identity(), Vec6::Zero());
	state0.g.block(0,3,3,1)=(Vec3()<<0.1,0.1,0.1).finished();
	quadrotor::State state_ref(Mat4::Identity(), Vec6::Zero());
	std::list<quadrotor::State> list_ref(2000, state_ref);
	std::list<quadrotor::U> list_u0(2000,Vec4::Ones());

	Mat12 M=Mat12::Identity();
	M.block(3,3,3,3)*=3;
	M.block(6,6,6,6)*=1e-1;
	Mat4 R=Mat4::Identity()*1e-2;
	Mat12 Mf=100*M;

	Vec4 umin=Vec4::Zero();
	Vec4 umax=Vec4::Ones()*6;

	quadrotor::iLQG(   sys,      dt,        T,
			        state0, list_u0, list_ref,
					     M,       R,       Mf,
					  umin,    umax);

/*************************************************************************
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

	quadrotor::State state=list_ref.front();

	std::list<quadrotor::U>::const_iterator itr_u=list_u0.begin();
	std::list<quadrotor::State> list_state;
	list_state.push_back(state);

	double dt=0.01;

	for(int n=0;n<1800;n++)
	{
		quadrotor::DState dstate(sys,state,*itr_u);
		state=state.update(dstate,dt);
		itr_u++;
		list_state.push_back(state);
	}
*************************************************************************/

//	std::string path="~/Documents/C++/test/bin/states.mat";
//	quadrotor::State::save(list_state,path);

}
