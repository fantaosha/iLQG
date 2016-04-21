/*************************************************************************
    > File Name: test.cpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Wed 20 Apr 2016 04:38:23 PM CDT
 ************************************************************************/
#include <Eigen/Dense>
#include <snOPT.hpp>

#include <type.hpp>
int main()
{
	Mat4 H1=(Vec4()<<3.2,1.1,5.3,2.9).finished().asDiagonal();
	Vec4 h1=(Vec4()<<-2.2,1.9,4.7,-3.1).finished();
	Vec4 x1=Vec4::Ones();
	Vec4 xmul1=Vec4::Zero();
	Eigen::Matrix<int,4,1> xstate1=Eigen::Matrix<int,4,1>::Zero();
	Vec4 xmin1=Vec4::Zero();
	Vec4 xmax1=Vec4::Ones()*2;
	double F1=0, Fmul1=0;
	int Fstate1=0;
	snOPT<4>::Params params1(xmin1,xmax1);



	Mat6 H2=(Vec6()<<1.7, 3.1, 6.7, 2.2, 9.1, 1.3).finished().asDiagonal();
	Vec6 h2=(Vec6()<<-1.7, 2.9, -3.7, 4.1, 8.9, -6.1).finished();
	Vec6 x2=Vec6::Ones();
	Vec6 xmul2=Vec6::Zero();
	Eigen::Matrix<int,6,1> xstate2=Eigen::Matrix<int,6,1>::Zero();
	Vec6 xmin2=Vec6::Zero();
	Vec6 xmax2=Vec6::Ones()*2;
	double F2=0, Fmul2=0;
	int Fstate2=0;
	snOPT<6>::Params params2(xmin2,xmax2);
	
	snOPT<4>::init();
	snOPT<4>::fmin(H1,    h1,  params1,
				   x1, xmul1,  xstate1,
			       F1, Fmul1,  Fstate1);

	snOPT<6>::init();
	snOPT<6>::fmin(H2,    h2,  params2,
				   x2, xmul2,  xstate2,
			       F2, Fmul2,  Fstate2);

	snOPT<4>::fmin(H1,    h1,  params1,
				   x1, xmul1,  xstate1,
			       F1, Fmul1,  Fstate1);
}
