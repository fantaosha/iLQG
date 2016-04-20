/*************************************************************************
    > File Name: cost.hpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Sun 17 Apr 2016 12:37:08 AM CDT
 ************************************************************************/
#ifndef QUADCOST
#define QUADCOST
#include <Eigen/Dense>
#include <liegroup.hpp>

namespace quadcost 
{

	double L(const Mat12 &M, const Mat4 &R, const Vec12 &dg, const Vec4 &du)
	{

		return (dg.transpose()*M*dg+du.transpose()*R*du)(0)*0.5;
	}

	Vec12 Lx(const Mat12 &M, const Vec12 &dg)
	{
		Vec12 Lx=M*dg;
		Lx.head(6)=SE3::dexpinv(-dg.head(6))*Lx.head(6);

		return Lx;
	}

	Mat12 Lxx(const Mat12 &M, const Vec12 &dg)
	{
		Vec6 dg_x=dg.head(6);
		Vec6 dg_v=dg.tail(6);

		Mat6 dexpinv=SE3::dexpinv(-dg_x);
		Mat6 dexpinvT=dexpinv.transpose();

		Eigen::Matrix<double,6,36> ddexpinv=SE3::ddexpinv(-dg_x);

		Mat12 Lxx=M;
		Lxx.block(0,0,6,6)=dexpinvT*M.block(0,0,6,6)*dexpinv;
		Lxx.block(0,6,6,6)=dexpinvT*M.block(0,6,6,6);
		Lxx.block(6,0,6,6)=M.block(0,6,6,6).transpose();

		Eigen::Matrix<double,1,6> r1=dg.transpose()*M.block(0,0,12,6);

		Mat6 DM1=-dexpinvT*(Mat6()<<r1*ddexpinv.block(0,0,6,6),
									r1*ddexpinv.block(0,6,6,6),
									r1*ddexpinv.block(0,12,6,6),
									r1*ddexpinv.block(0,18,6,6),
									r1*ddexpinv.block(0,24,6,6),
									r1*ddexpinv.block(0,30,6,6)).finished();
		
		Eigen::Matrix<double,1,6> r2=-0.5*r1*dexpinv;

		Mat6 DM2=(Mat6()<<r2*SE3::ad(SE3::e[0]),
						  r2*SE3::ad(SE3::e[1]),
						  r2*SE3::ad(SE3::e[2]),
						  r2*SE3::ad(SE3::e[3]),
						  r2*SE3::ad(SE3::e[4]),
						  r2*SE3::ad(SE3::e[5])).finished();
		
		Lxx.block(0,0,6,6)+=DM1+DM2;	

		return Lxx;
	}

	Vec4 Lu(const Mat4 &R, const Vec4 &du)
	{
		return R*du;
	}

	Mat4 Luu(const Mat4 &R, const Vec4 &du)
	{
		return R;
	}
}
#endif
