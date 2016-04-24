/*************************************************************************
    > File Name: car_ddp.cpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Sat 23 Apr 2016 07:48:48 PM CDT
 ************************************************************************/
#include <car.hpp>
#include <Eigen/Dense>

#include <car.hpp>
#include <ddp.hpp>

int main()
{
	car::System sys;
	double dt=0.01;
	double N=800;

	Vec5 x0=(Vec5()<<3,3,4.5,0,0).finished();
	car::State state0(x0);
	car::State xref(Vec5::Zero());

	std::vector<car::U> us0(2000,Vec2::Zero());
	std::vector<car::State> xrefs(2000, xref);

	Mat5 M=(Vec5()<<0.3,0.3,0.3,0.1,0.1).finished().asDiagonal();
	Mat5 Mf=M*50;
	Mat2 R=(Vec2()<<0.1,0.1).finished().asDiagonal();

	DDP<car>::Params params(M,R,Mf);

	Vec2 umin;
	umin(0)=-0.5;
	umin(1)=-2;
	Vec2 umax;
	umax(0)=0.5;
	umax(1)=2;

	DDP<car> ddp(sys,dt);
	std::vector<Vec2> us;

	timespec T_start, T_end;
	clock_gettime(CLOCK_MONOTONIC,&T_start);
	ddp.init(state0, us0, xrefs, params, umin, umax, N);
	ddp.iterate(10, us);
	clock_gettime(CLOCK_MONOTONIC,&T_end);
	std::cout<<"time consumed is "<<(T_end.tv_sec-T_start.tv_sec)+(T_end.tv_nsec-T_start.tv_nsec)/1000000000.0<<"s"<<std::endl;

	return 0;
}
