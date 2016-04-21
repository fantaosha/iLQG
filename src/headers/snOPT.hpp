/*************************************************************************
    > File Name: qopt.hpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Mon 18 Apr 2016 02:03:16 PM CDT
 ************************************************************************/
#ifndef _SNOPT_
#define _SNOPT_
#include <stdio.h>
#include <string.h>
#include <iostream>

#include <Eigen/Dense>
#include <snoptProblem.hpp>

template<size_t N> class snOPT
{
	public:
		static const size_t num=N;
		static const size_t nef=1;

		typedef Eigen::Matrix<double,N,N> Mat;
		typedef Eigen::Matrix<double,N,1> Vec;

		static const int Cold=0;
		static const int Basis=1;
		static const int Warm=2;

		static const int lenG=N;
		static const int neG=N;
		static const int lenA=N;
		static const int neA=0;

		const static int ObjRow=0;

		constexpr static double inf=1.0e20;
		
		static Mat H;
		static Vec h;

		struct Params
		{
			public:
				Vec xlow;
				Vec xupp;
				double Flow;
				double Fupp;

				double ObjAdd;

				Params():xlow(-Vec::Ones()*1e20), xupp(Vec::Ones()*1e20), Flow(-1e20), Fupp(1e20), ObjAdd(0)
				{
				}

				Params(Vec xlow_, Vec xupp_, double Flow_=-1e20, double Fupp_=1e20, double ObjAdd_=0):xlow(xlow_), xupp(xupp_), Flow(Flow_), Fupp(Fupp_ ), ObjAdd(ObjAdd_)
				{
				}
		};

	protected:
		static snoptProblemA problem;

		static bool initialized;
			   
		static int iGfun[lenG];
		static int jGvar[lenG];

		static int iAfun[lenA];
		static int jAvar[lenA];
		static double A[lenA];

		static double x[num];
		static double xlow[num];
		static double xupp[num];
		static double xmul[num];
		static int xstate[num];

		static double F[nef];
		static double Flow[nef];
		static double Fupp[nef];
		static double Fmul[nef];
		static int Fstate[nef];

		static double ObjAdd;

	public:
		static void init()
		{
			if (!initialized)
			{
				problem.setProbName("Optimization");		
				problem.setProblemSize(num,nef);

				problem.setIntParameter( "Major print level", 0);
				problem.setIntParameter( "Print frequency", 0);
				problem.setIntParameter( "Summary file", 0);
				problem.setIntParameter( "Summary frequency", 0);
				problem.setIntParameter( "Hessian updates", 5);
				problem.setIntParameter( "Minor print level", 0 );
				problem.setIntParameter( "Derivative option", 3 );
				problem.setIntParameter( "Verify level ", 0 );
				problem.setIntParameter( "Major Iteration limit", 20 );
				problem.setIntParameter( "Iterations limit", 200 );
				problem.setRealParameter( "Major optimality tolerance", 5e-4 );

				for(size_t j=0;j<num;j++)
					jGvar[j]=j;

				problem.setG(lenG, neG, iGfun, jGvar);
				problem.setA(lenA, neA, iAfun, jAvar,A);
				problem.setUserFun(usrfunG);

				initialized=true;
			}
		}
	
		static int is_ready()
		{
			return initialized;
		}

		void static usrfunObj(int *Status,  int     *n,  double  x[],
					int  *needF,  int   *neF,  double  F[],
					int  *needG,  int   *neG,  double  G[],
					char    *cu,  int *lencu,  int    iu[], 
					int  *leniu, double ru[],  int  *lenru)
		{
			Vec xx;
			memcpy(xx.data(),x,num*sizeof(double));
			Eigen::Matrix<double,1,num> xxT=xx.transpose();

			F[0]=(0.5*xxT*H*xx+xxT*h)(0);
		}

		void static usrfunG(int *Status,  int     *n,  double  x[],
					int  *needF,  int   *neF,  double  F[],
					int  *needG,  int   *neG,  double  G[],
					char    *cu,  int *lencu,  int    iu[], 
					int  *leniu, double ru[],  int  *lenru)
		{
			Vec xx;
			memcpy(xx.data(),x,num*sizeof(double));
			Eigen::Matrix<double,1,num> xxT=xx.transpose();
			if(*needF>0)
				F[0]=(0.5*xxT*H*xx+xxT*h)(0);
			
			if(*needG>0)
			{
				Vec dx=H*xx+h;
				memcpy(G, dx.data(), sizeof(double)*num);
			}

		}

		static int fmin(Mat const &  H, Vec const & h, Params const & params, 
				  Vec & x0, Vec & xmul0, Eigen::Matrix<int,N,1> & xstate0,
				  double & F0, double & Fmul0, int & Fstate0)
		{
			snOPT::H=H;
			snOPT::h=h;

			memcpy(x,x0.data(),sizeof(double)*num);
			memcpy(xmul,xmul0.data(),sizeof(double)*num);
			memcpy(xstate,xstate0.data(),sizeof(int)*num);
			
			F[0]=F0;
			Fmul[0]=Fmul0;
			Fstate[0]=Fstate0;
		
			memcpy(xlow,params.xlow.data(), sizeof(double)*num);
			memcpy(xupp,params.xupp.data(), sizeof(double)*num);
			Flow[0]=params.Flow;
			Fupp[0]=params.Fupp;

			ObjAdd=params.ObjAdd;


			problem.setObjective  ( ObjRow, ObjAdd );
			problem.setX          ( x, xlow, xupp, xmul, xstate );
			problem.setF          ( F, Flow, Fupp, Fmul, Fstate );

			int result=problem.solve          ( Warm );

			memcpy(x0.data(),x,sizeof(double)*num);
			memcpy(xmul0.data(),xmul,sizeof(double)*num);
			memcpy(xstate0.data(),xstate,sizeof(int)*num);

			F0=F[0];
			Fmul0=Fmul[0];
			Fstate0=Fstate[0];
			return result;
		}
};

template <size_t N> int snOPT<N>::iGfun[]={0};
template <size_t N> int snOPT<N>::jGvar[]={0};
template <size_t N> int snOPT<N>::iAfun[]={0};
template <size_t N> int snOPT<N>::jAvar[]={0};
template <size_t N> double snOPT<N>::A[]={0};

template <size_t N> double snOPT<N>::x[snOPT<N>::num]={0};
template <size_t N> double snOPT<N>::xlow[snOPT<N>::num]={0};
template <size_t N> double snOPT<N>::xupp[snOPT<N>::num]={0};
template <size_t N> double snOPT<N>::xmul[snOPT<N>::num]={0};
template <size_t N> int snOPT<N>::xstate[snOPT<N>::num]={0};

template <size_t N> double snOPT<N>::F[snOPT<N>::nef]={0};
template <size_t N> double snOPT<N>::Flow[snOPT<N>::nef]={0};
template <size_t N> double snOPT<N>::Fupp[snOPT<N>::nef]={0};
template <size_t N> double snOPT<N>::Fmul[snOPT<N>::nef]={0};
template <size_t N> int snOPT<N>::Fstate[snOPT<N>::nef]={0};

template <size_t N> double snOPT<N>::ObjAdd=0;
template <size_t N> bool snOPT<N>::initialized=false;

template <size_t N> typename snOPT<N>::Mat snOPT<N>::H=snOPT<N>::Mat::Identity();
template <size_t N> typename snOPT<N>::Vec snOPT<N>::h=snOPT<N>::Vec::Zero();
template <size_t N> snoptProblemA snOPT<N>::problem;
//Vec quadopt::h=Vec4::Zero();

#endif
