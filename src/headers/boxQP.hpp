/*************************************************************************
    > File Name: boxQP.hpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Fri 22 Apr 2016 03:20:23 PM CDT
 ************************************************************************/
#ifndef _BOXQP
#define _BOXQP
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>

template<size_t N> Eigen::Matrix<double,N,1> clamp(Eigen::Matrix<double,N,1> const & x, 
		                      Eigen::Matrix<double,N,1> const & xmin, Eigen::Matrix<double,N,1> const & xmax)
{
	Eigen::Matrix<double,N,1> xn=x;
	for(int i=0;i<N;i++)
		if(xn(i)<xmin(i))
			xn(i)=xmin(i);
		else
			if(xn(i)>xmax(i))
				xn(i)=xmax(i);

	return xn;
}

template<size_t N> Eigen::MatrixXd select(Eigen::Matrix<double,N,N> const & M, Eigen::Array<int,N,1> const & istate, Eigen::Array<int,N,1> const & jstate)
{
	
	int m=istate.sum();
	int n=jstate.sum();

	Eigen::MatrixXd M1=Eigen::MatrixXd::Zero(m,N);
	Eigen::MatrixXd M2=Eigen::MatrixXd::Zero(m,n);
	
	for(int i=0,j=0;i<N && j<m;i++)
		if(istate(i))
			M1.row(j++)=M.row(i);

	for(int i=0,j=0; i<N && j<n; i++)
		if(jstate(i))
			M2.col(j++)=M1.col(i);

	return M2;
}

template<size_t N> Eigen::VectorXd select(Eigen::Matrix<double,N,1> const & v, Eigen::Array<int,N,1> const & istate)
{
	
	int n=istate.sum();

	Eigen::MatrixXd v1=Eigen::MatrixXd::Zero(n,1);
	
	for(int i=0,j=0;i<N && j<n;i++)
		if(istate(i))
			v1(j++)=v(i);

	return v1;
}

template<size_t N> void assign(Eigen::Matrix<double,N,1> &v1, Eigen::MatrixXd const & v2, Eigen::Array<int,N,1> const & istate)
{
	assert(v2.rows()>=istate.sum());
	for(int i=0,j=0;i<N;i++)
		if(istate(i))
			v1(i)=v2(j++);
}

template<size_t N> void boxQP(Eigen::Matrix<double,N,N> const & H,    Eigen::Matrix<double,N,1> const & g, 
		                      Eigen::Matrix<double,N,1> const & xmin, Eigen::Matrix<double,N,1> const & xmax, 
							  Eigen::Matrix<double,N,1> const & x0,
							  Eigen::Matrix<double,N,1>       &  x, 
							  Eigen::MatrixXd                 & Hfree, Eigen::Array<int,N,1> & free,
							  int & result)
{
	typedef Eigen::Matrix<double,N,N> MatNd;
	typedef Eigen::Matrix<double,N,1> VecNd;
	typedef Eigen::Array<int,N,1> IndexN;
	typedef Eigen::MatrixXd Mat;
	typedef Eigen::VectorXd Vec;

	// Define control parameters
	constexpr static size_t maxIter=100;
	constexpr static double minGrad=1e-8;
	constexpr static double relTol=1e-8;
	constexpr static double stepDec=0.6;
	constexpr static double minStep=1e-22;
	constexpr static double Armijo=0.1;

	IndexN clamped=IndexN::Zero();
	free=1-clamped;

#ifdef DEBUG
	std::cout<<"============================================================="<<std::endl;
	std::cout<<"H: "<<std::endl<<H<<std::endl;
	std::cout<<"g: "<<g.transpose()<<std::endl;
	std::cout<<"xmin: "<<xmin.transpose()<<std::endl;
	std::cout<<"xmax: "<<xmax.transpose()<<std::endl;
	std::cout<<"x0: "<<x0.transpose()<<std::endl;
	std::cout<<"============================================================="<<std::endl;
#endif
	x=clamp<N>(x0,xmin,xmax);

	result=0;
	double gnorm;

	double J=x.dot(0.5*H*x+g);
	double J0=J;
	
	double gnom;

	size_t iter;

	Eigen::LLT<Mat> llt;

	for(iter=0; iter<maxIter; iter++)
	{
		if(result)
			break;

		if(iter>0 && fabs(J0-J)<relTol*J0)
		{
			result=4;
			break;
		}
		J0=J;
		VecNd grad=H*x+g;
		IndexN clamped0=clamped;
		clamped=IndexN::Zero();
		
		for(int i=0;i<N;i++)	
			clamped(i)= (x(i)<=xmin(i) && grad(i)>0) || (x(i)>=xmax(i) && grad(i)<0);
		
		free=1-clamped;

		if(clamped.sum()==N)
		{
			result=6;
			break;
		}
		
		if(iter==0 || (clamped-clamped0).abs().sum()>0)
		{
			Mat Hs = select<N>(H,free,free);	
			llt.compute(Hs);

			if(llt.info()!=Eigen::Success)
			{
				result=-1;
				break;
			}

			Hfree=llt.matrixU();
		}

		gnorm=select<N>(grad,free).norm();

		if(gnorm<minGrad)
		{
			result=5;
			break;
		}

		VecNd grad_clamped=g+H*((x.array()*(clamped.template cast<double>())).matrix());
		VecNd dx=VecNd::Zero();
		assign<N>(dx,-llt.solve(select<N>(grad_clamped,free))-select<N>(x,free),free);
		
#ifdef DEBUG
		std::cout<<"========================================="<<std::endl;
		std::cout<<"Iteration #"<<iter+1<<std::endl;
		std::cout<<"========================================="<<std::endl;
		std::cout<<"x0:"<<x.transpose()<<std::endl;
		std::cout<<"J0: "<<J0<<std::endl;
		std::cout<<"gradient: "<<grad.transpose()<<std::endl;
		std::cout<<"clamped index: "<<clamped.transpose()<<std::endl;
		std::cout<<"gradient clamped: "<<grad_clamped.transpose()<<std::endl;
		std::cout<<"dx: "<<dx.transpose()<<std::endl;
#endif

		double dxdotg=dx.dot(grad);
		if(dxdotg>=0)
			break;

		double step=1; 
		VecNd xc=clamp<N>(x+step*dx, xmin, xmax);
		double Jn=xc.dot(0.5*H*xc+g);

		while((Jn-J0)/(step*dxdotg)<Armijo)
		{
			
			step*=stepDec;
			xc=clamp<N>(x+step*dx, xmin, xmax);
			Jn=xc.dot(0.5*H*xc+g);
		
			if(step<minStep)
			{
				result=2;
				break;
			}
		}

		x=xc;
		J=Jn;

#ifdef DEBUG
		std::cout<<"xn "<<x.transpose()<<std::endl;
		std::cout<<"Jn "<<Jn<<std::endl;
#endif
	}

	if(iter>=maxIter) 
		result=1;
	
}

#endif
