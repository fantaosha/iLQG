/*************************************************************************
    > File Name: qopt.hpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Mon 18 Apr 2016 02:03:16 PM CDT
 ************************************************************************/
#include <stdio.h>
#include <string.h>
#include <iostream>

#include <Eigen/Dense>
#include <snoptProblem.hpp>
#include <type.hpp>

class quadopt
{
	public:
		const static size_t num=4;
		const static size_t nef=1;
		
		const static int Cold=0;
		const static int Basis=1;
		const static int Warm=2;

		const static int lenG=4;
		const static int neG=4;
		const static int lenA=4;
		const static int neA=0;

		const static int ObjRow=0;

		constexpr static double inf=1.0e20;


		static Mat4 H;
		static Vec4 h;

		struct Params
		{
			public:
				Vec4 xlow;
				Vec4 xupp;
				double Flow;
				double Fupp;

				double ObjAdd;

				Params():xlow(-Vec4::Ones()*1e20), xupp(Vec4::Ones()*1e20), Flow(-1e20), Fupp(1e20), ObjAdd(0)
				{
				}

				Params(Vec4 xlow_, Vec4 xupp_, double Flow_=-1e20, double Fupp_=1e20, double ObjAdd_=0):xlow(xlow_), xupp(xupp_), Flow(Flow_), Fupp(Fupp_ ), ObjAdd(ObjAdd_)
				{
				}
		};

	protected:
		static snoptProblemA problem;
			   
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
		quadopt()
		{
			problem.setProbName("Quadratic Optimization");		
//			problem.setPrintFile("quadrotor.out");		
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

			problem.setG(lenG, neG, iGfun, jGvar);
			problem.setA(lenA, neA, iAfun, jAvar,A);
			problem.setUserFun(usrfunG);
		}
		
		void static usrfunObj(int *Status,  int     *n,  double  x[],
					int  *needF,  int   *neF,  double  F[],
					int  *needG,  int   *neG,  double  G[],
					char    *cu,  int *lencu,  int    iu[], 
					int  *leniu, double ru[],  int  *lenru)
		{
			Vec4 xx;
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
			Vec4 xx;
			memcpy(xx.data(),x,num*sizeof(double));
			Eigen::Matrix<double,1,num> xxT=xx.transpose();
			if(*needF>0)
				F[0]=(0.5*xxT*H*xx+xxT*h)(0);
			
			if(*needG>0)
			{
				Vec4 dx=H*xx+h;
				memcpy(G, dx.data(), sizeof(double)*num);
			}

		}

		int fmin(Mat4 const &  H, Vec4 const & h, Params const & params, 
				  Vec4 & x0, Vec4 & xmul0, Eigen::Matrix<int,4,1> & xstate0,
				  double & F0, double & Fmul0, int & Fstate0)
		{
			this->H=H;
			this->h=h;

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

int quadopt::iGfun[]={0,0,0,0};
int quadopt::jGvar[]={0,1,2,3};

int quadopt::iAfun[]={0,0,0,0};
int quadopt::jAvar[]={0,1,2,3};
double quadopt::A[]={0};

snoptProblemA quadopt::problem;

double quadopt::x[num]={0};
double quadopt::xlow[num]={0};
double quadopt::xupp[num]={0};
double quadopt::xmul[num]={0};
int quadopt::xstate[num]={0};

double quadopt::F[nef]={0};
double quadopt::Flow[nef]={0};
double quadopt::Fupp[nef]={0};
double quadopt::Fmul[nef]={0};
int quadopt::Fstate[nef]={0};

double quadopt::ObjAdd=0;

Mat4 quadopt::H=Mat4::Identity();
Vec4 quadopt::h=Vec4::Zero();
