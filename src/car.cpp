/*************************************************************************
    > File Name: car.cpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Thu 21 Apr 2016 06:55:51 AM CDT
 ************************************************************************/
#include <car.hpp>
#include <Eigen/Dense>

#include <car.hpp>
#include <ilqg.hpp>

int main()
{
	car::System sys;
	double dt=0.1;
	double N=33;

	Vec5 x0=(Vec5()<<-5,-2,-1.2,0,0).finished();
	car::State state0(x0);
	car::State state_ref(Vec5::Zero());

	std::list<car::U> list_u0(2000,Vec2::Zero());
	std::list<car::State> list_ref(2000, state_ref);

	Mat5 M=(Vec5()<<5,5,1,1,1).finished().asDiagonal();
	Mat5 Mf=M;
	Mat2 R=(Vec2()<<1,5).finished().asDiagonal();

	iLQG<car>::Params params(M,R,Mf);

	Vec2 umin=Vec2::Ones()*4;
	Vec2 umax=Vec2::Ones()*4;
	
	iLQG<car> ilqg(sys,dt);
	ilqg.init(state0, list_u0, list_ref, params, umin, umax,N);
	ilqg.evaluate(-1,20,list_u0);
	return 0;
}
