/*************************************************************************
    > File Name: liegroup.hpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Fri 20 Feb 2015 07:15:26 PM EST
 ************************************************************************/
#ifndef LIEGROUP
#define LIEGROUP

#include<Eigen/Dense>
#include<type.hpp>
#include<cmath>

//#include<iostream>

#ifndef PI
#define PI 3.14159265358979323846
#endif

class SO2
{
	public:
		static const uchar row;
		static const uchar col;
		static const uchar dim;

		static const Mat2 E[1]; 

		static Mat2 exp(double th);
		static double log(Mat2 R);

		static Mat2 hat(double th);
		static double vee(Mat2);
};

const uchar SO2::row=2;
const uchar SO2::col=2;
const uchar SO2::dim=1;

const Mat2 SO2::E[]={(Mat2()<<0,-1,1,0).finished()};

Mat2 SO2::exp(double th)
{
	return (Mat2()<<cos(th), -sin(th), sin(th), cos(th)).finished();
}

double SO2::log(Mat2 R)
{
	return atan2(R(1,0), R(0,0));
}

Mat2 SO2::hat(double th)
{
	return (Mat2()<<0,-th,th,0).finished();
}

double SO2::vee(Mat2 W)
{
	return W(1,0);
}

class SE2
{
	public:
		static const uchar row;
		static const uchar col;
		static const uchar dim;

		static const Vec3 e[3]; 
		static const Mat3 E[3]; 


		static Mat3 exp(Vec3 v);
		static Vec3 log(Mat3 g);
		static Mat3 dexp(Vec3 v);
		static Mat3 dexpinv(Vec3 v);
		static Mat3 dexpinvAd(Vec3 v);

		static Mat3 ad(Vec3 v);
		static Mat3 Ad(Mat3 g);

		static Mat3 hat(Vec3 v);
		static Vec3 vee(Mat3);

		static Mat3 cay(Vec3 v);
		static Vec3 cayinv(Mat3 g);
		static Mat3 dcay(Vec3 v);
		static Mat3 dcayinv(Vec3 v);

		static Mat3 inverse(Mat3 g);

};

const uchar SE2::row=3;
const uchar SE2::col=3;
const uchar SE2::dim=3;


const Vec3 SE2::e[]={(Vec3()<<1,0,0).finished(),
					 (Vec3()<<0,1,0).finished(), 
					 (Vec3()<<0,0,1).finished()};

const Mat3 SE2::E[]={(Mat3()<<0,-1,0,1,0,0,0,0,0).finished(),
					 (Mat3()<<0,0,1,0,0,0,0,0,0).finished(), 
					 (Mat3()<<0,0,0,0,0,1,0,0,0).finished()};

Mat3 SE2::exp(Vec3 v)
{
	Mat3 g=Mat3::Identity();

	if(fabs(v(0))<1e-12)
	{
		g(0,2)=v(1);
		g(1,2)=v(2);
	}
	else
	{
		Mat2 R((Mat2()<<cos(v(0)), -sin(v(0)), sin(v(0)), cos(v(0))).finished());
		Vec2 p((Vec2()<<v(2),-v(1)).finished());
		g.block(0,0,2,2)=R;
		g.block(0,2,2,1)=((R-Mat2::Identity())/v(0))*p;
	}

	return g;
}

Vec3 SE2::log(Mat3 g)
{
	double th=atan2(g(1,0),g(0,0));
	Vec3 v;
	if(fabs(th)<4e-6)
		v<<th,g(0,2),g(1,2);
	else
	{
		double a=g(1,0)*(th/(g(0,0)-1));
		v<<th,((Mat2()<<-a, th, -th, -a).finished()*0.5)*g.block(0,2,2,1);
	}

	return v;
}

Mat3 SE2::dexp(Vec3 v)
{
	if(fabs(v(0))<1e-6)
		return Mat3::Identity()+ad(v)/2;
	else
	{
		Mat3 dexp=Mat3::Zero();

		double sth=sin(v(0));
		double cth=cos(v(0));

		dexp<<                                         1,            0,             0,
			  (v(1)*(v(0)-sth)+v(2)*(1-cth))/(v(0)*v(0)),     sth/v(0), -(1-cth)/v(0),
			  (v(2)*(v(0)-sth)-v(1)*(1-cth))/(v(0)*v(0)), (1-cth)/v(0),      sth/v(0);
		
		return dexp;
	}
}

Mat3 SE2::dexpinv(Vec3 v)
{
	if(fabs(v(0))<1e-6)
		return Mat3::Identity()-ad(v)/2;
	else
	{
		Mat3 dexpinv=Mat3::Zero();

		double sth=sin(v(0));
		double cth=cos(v(0));

		dexpinv<<                                               1,                     0,                    0,
			  -v(2)/2+v(1)*(v(0)*sth+2*cth-2)/(2*v(0)*(cth-1)),  v(0)*sth/(2*(1-cth)),                v(0)/2,
			   v(1)/2+v(2)*(v(0)*sth+2*cth-2)/(2*v(0)*(cth-1)),               -v(0)/2,  v(0)*sth/(2*(1-cth));

		return dexpinv;
	}
}

Mat3 SE2::dexpinvAd(Vec3 v)
{
//	return dexpinv(v)*Ad(exp(v));
	return dexpinv(-v);
}

Mat3 SE2::ad(Vec3 v)
{
	return (Mat3()<<   0,    0,     0,
					v(2),    0, -v(0),
				   -v(1), v(0),   0).finished();
}

Mat3 SE2::Ad(Mat3 g)
{
	Mat3 gAd=Mat3::Identity();
	gAd.block(1,1,2,2)=g.block(0,0,2,2);
	gAd.block(1,0,2,1)<<g(1,2),-g(0,2);

	return gAd;
}

Mat3 SE2::hat(Vec3 v)
{
	return (Mat3()<<0, -v(0), v(1),
				    v(0),  0, v(2),
					0,     0,    0).finished();
}

Vec3 SE2::vee(Mat3 V)
{
	return (Vec3()<<V(1,0), V(0,2), V(1,2)).finished();
}

Mat3 SE2::cay(Vec3 v)
{
	double w2=v(0)*v(0);
	
	Mat3 g=Mat3::Identity();
	g.block(0,0,2,3)<<4-w2,  -4*(v(0)),  -2*(v(0)*v(2)-2*v(1)),
					  4*v(0),   4-w2,   2*(v(0)*v(1)+2*v(2));
	g.block(0,0,2,3)/=4+w2;

	return g;
		
}

Vec3 SE2::cayinv(Mat3 g)
{
	Vec3 v;
	v(0)=2.0*g(1,0)/(1+g(0,0));
	v.tail(2)=(Mat2()<<1, v(0)/2, -v(0)/2, 1).finished()*g.block(0,2,2,1);

	return v;
}

Mat3 SE2::dcay(Vec3 v)
{
	return 4.0/(4+v(0)*v(0))*(Mat3::Identity()+0.5*SE2::ad(v));
}

Mat3 SE2::dcayinv(Vec3 v)
{
	return Mat3::Identity()-SE2::ad(v)/2.0+ (Mat3()<<v(0)*v, Eigen::Matrix<double,3,2>::Zero()).finished()/4.0;
}


Mat3 SE2::inverse(Mat3 g)
{
	Mat3 gi=Mat3::Identity();
	gi.block(0,0,2,2)=g.block(0,0,2,2).transpose();
	gi.block(0,2,2,1)=-gi.block(0,0,2,2)*g.block(0,2,2,1);

	return gi;

}

class SO3
{
	public:
		static const uchar row;
		static const uchar col;
		static const uchar dim;

		static const Vec3 e[3];
		static const Mat3 E[3];


		static Mat3 exp(Vec3 v);
		static Vec3 log(Mat3 g);
		static Mat3 dexp(Vec3 v);
		static Mat3 dexpinv(Vec3 v);
		static Mat3 dexpinvAd(Vec3 v);
		static Eigen::Matrix<double,3,9> ddexp(Vec3);
		static Eigen::Matrix<double,3,9> ddexpinv(Vec3);


		static Mat3 cay(Vec3 v);
		static Vec3 cayinv(Mat3 g);
		static Mat3 dcay(Vec3 v);
		static Mat3 dcayinv(Vec3 v);
		static Mat3 dcayinvAd(Vec3 v);
		static Eigen::Matrix<double,3,9> ddcayinv(Vec3 v, double dt);
		static Eigen::Matrix<double,3,9> ddcayinvAd(Vec3 v, double dt);

		static Mat3 hat(Vec3 v);
		static Vec3 vee(Mat3 vh);

		static Mat3 Ad(Mat3 g);	
		static Mat3 ad(Vec3 v);

		static Mat3 inverse(Mat3 g);

		static Mat4 toSE3(Mat3 g);
		static Vec6 tose3(Vec3 v);

		static Vec3 getRPY(Mat3);
};

const uchar SO3::row=3;
const uchar SO3::col=3;
const uchar SO3::dim=3;

const Vec3 SO3::e[]={(Vec3()<<1.0,0.0,0.0).finished(), 
					 (Vec3()<<0.0,1.0,0.0).finished(), 
					 (Vec3()<<0.0,0.0,1.0).finished()};

const Mat3 SO3::E[]={(Mat3()<<0.0,0.0,0.0,0.0,0.0,-1,0.0,1,0.0).finished(), 
					 (Mat3()<<0.0,0.0,1,0.0,0.0,0.0,-1,0.0,0.0).finished(), 
					 (Mat3()<<0.0,-1,0.0,1,0.0,0.0,0.0,0.0,0.0).finished()};

Mat3 SO3::Ad(Mat3 g)
{
	return g;
}

Mat3 SO3::ad(Vec3 V)
{
	return hat(V);
}


Mat3 SO3::exp(Vec3 V)
{
	Mat3 Vh=hat(V);
	double th=V.norm();
	if (th<1e-10)
		return Mat3::Identity();
	else
	{
		return Mat3::Identity()+sin(th)*Vh/th+(1-cos(th))*Vh*Vh/(th*th);
	}
}

Vec3 SO3::log(Mat3 g)
{
	Vec3 w=Vec3::Zero();
	double TR=g.trace();
	if(fabs(TR-3)>1e-10) 
		if(fabs(TR+1)<5e-10)
		{
			double w11=(1+g(0,0))/2.0;
			double w22=(1+g(1,1))/2.0;
			double w33=(1+g(2,2))/2.0;

			if(w11>w22)
				if(w11>w33)
				{
					double w1=sqrt(w11);
					w<<w1,g(0,1)/(2*w1),g(0,2)/(2*w1);
				}
				else
				{
					double w3=sqrt(w33);
					w<<g(2,0)/(2*w3),g(2,1)/(2*w3),w3;
				}
			else
				if(w22>w33)
				{
					double w2=sqrt(w22);
					w<<g(1,0)/(2*w2),w2,g(1,2)/(2*w2);
				}
				else
				{
					double w3=sqrt(w33);
					w<<g(2,0)/(2*w3),g(2,1)/(2*w3),w3;
				}

			w*=PI;
		}
		else
		{
			double th=acos((TR-1)/2);
			w(0)=g(2,1)-g(1,2);
			w(1)=g(0,2)-g(2,0);
			w(2)=g(1,0)-g(0,1);
			w=w*th/(2*sin(th));
		}

	return w;
}

Mat3 SO3::dexp(Vec3 V)
{
	Mat3 Vh=hat(V);
	double th=V.norm();

	if (th<1e-6)
		return Mat3::Identity();
	else
	{
		double sth=sin(th);
		double cth=cos(th);

		double th2=th*th;
		double th3=th2*th;

		return Mat3::Identity()+((1-cth)/th2)*Vh+((th-sth)/th3)*Vh*Vh;
	}
}

Mat3 SO3::dexpinv(Vec3 V)
{
	Mat3 Vh=hat(V);
	double th=V.norm();

	if (th<1e-6)
		return Mat3::Identity();
	else
	{
		double sth=sin(th);
		double cth=cos(th);

		double th2=th*th;

		return Mat3::Identity()-Vh/2+((th*sth+2*cos(th)-2)/(2*th2*(cos(th)-1)))*Vh*Vh;
	}
}

Mat3 SO3::dexpinvAd(Vec3 V)
{
//	return dexpinv(V)*exp(V);
	return dexpinv(-V);
}

Eigen::Matrix<double,3,9> SO3::ddexp(Vec3 V)
{
	double th=V.norm();
	double th2=th*th;
	double th3=th2*th;
	double th4=th3*th;
	double th5=th4*th;

	Mat3 Vh=hat(V);

	double c1=0;
	double c2=0;

	double dc1=0;
	double dc2=0;

	double sth=sin(th);
	double cth=cos(th);

	Eigen::Matrix<double,3,9> DF;

	if (th>1e-6)
	{
		c1=(1-cth)/th2;
		c2=(th-sth)/th3;

		dc1=(2*(cth-1)+th*sth)/th4;
		dc2=-(2*th-3*sth+th*cth)/th5;

		Mat3 Dth=dc1*Vh+dc2*Vh*Vh;

		DF<<c1*E[0]+c2*(Vh*E[0]+E[0]*Vh)+Dth*V(0),
			c1*E[1]+c2*(Vh*E[1]+E[1]*Vh)+Dth*V(1),
			c1*E[2]+c2*(Vh*E[2]+E[2]*Vh)+Dth*V(2);
	}
	else
	{
		c1=1.0/2.0;
		c2=1.0/6.0;

		DF<<c1*E[0]+c2*(Vh*E[0]+E[0]*Vh),
			c1*E[1]+c2*(Vh*E[1]+E[1]*Vh),
			c1*E[2]+c2*(Vh*E[2]+E[2]*Vh);
	}

	return DF;
}

Eigen::Matrix<double,3,9> SO3::ddexpinv(Vec3 V)
{
	double th=V.norm();
	double th2=th*th;
	double th3=th2*th;
	double th4=th3*th;

	Mat3 Vh=hat(V);

	double c1=-0.5;
	double c2=0;

	double dc1=0;
	double dc2=0;
	
	double sth=sin(th);
	double cth=cos(th)-1;

	Eigen::Matrix<double,3,9> DF;

	if(th>=1e-6)
	{
		c2=(th*sth+2*cth)/(2*th2*cth);
		dc2=-(4*cth+th*sth+th2)/(2*th4*cth);

		Mat3 Dth=dc1*Vh+dc2*Vh*Vh;

		DF<<c1*E[0]+c2*(Vh*E[0]+E[0]*Vh)+Dth*V(0),
			c1*E[1]+c2*(Vh*E[1]+E[1]*Vh)+Dth*V(1),
			c1*E[2]+c2*(Vh*E[2]+E[2]*Vh)+Dth*V(2);
	}
	else
	{
		c2=-1.0/12.0;
		dc2=0;

		DF<<c1*E[0]+c2*(Vh*E[0]+E[0]*Vh),
			c1*E[1]+c2*(Vh*E[1]+E[1]*Vh),
			c1*E[2]+c2*(Vh*E[2]+E[2]*Vh);
	}


	return DF;
}

Mat3 SO3::cay(Vec3 V)
{
	Mat3 Vh=hat(V);
	double th=V.squaredNorm();
	return Mat3::Identity()+4/(4+th)*(Vh+Vh*Vh/2);
}

Vec3 SO3::cayinv(Mat3 g)
{
	Vec3 w=Vec3::Zero();
	w(0)=g(2,1)-g(1,2);
	w(1)=g(0,2)-g(2,0);
	w(2)=g(1,0)-g(0,1);
	
	return 2*w/(1.0+g.trace());
}

Mat3 SO3::dcay(Vec3 V)
{
	Mat3 Vh=hat(V);
	double th=V.squaredNorm();
	return 4/(4+th)*(Vh/2+Mat3::Identity());
}

Mat3 SO3::dcayinv(Vec3 V)
{
	Mat3 Vh=hat(V);
	return Mat3::Identity()-Vh/2+V*V.transpose()/4;
}

Mat3 SO3::dcayinvAd(Vec3 V)
{
	Mat3 Vh=hat(V);
	return Mat3::Identity()+Vh/2+V*V.transpose()/4;
}

Eigen::Matrix<double,3,9> SO3::ddcayinv(Vec3 V,double dt)
{
	Mat DF=Mat::Zero(3,9);
	DF<<-E[0]*dt/2+dt*dt*(e[0]*V.transpose()+V*e[0].transpose())/4,
		-E[1]*dt/2+dt*dt*(e[1]*V.transpose()+V*e[1].transpose())/4,
		-E[2]*dt/2+dt*dt*(e[2]*V.transpose()+V*e[2].transpose())/4;
	return DF;

}

Eigen::Matrix<double,3,9> SO3::ddcayinvAd(Vec3 V,double dt)
{
	Mat DF=Mat::Zero(3,9);
	DF<<E[0]*dt/2+dt*dt*(e[0]*V.transpose()+V*e[0].transpose())/4,
		E[1]*dt/2+dt*dt*(e[1]*V.transpose()+V*e[1].transpose())/4,
		E[2]*dt/2+dt*dt*(e[2]*V.transpose()+V*e[2].transpose())/4;
	return DF;
}

Mat3 SO3::hat(Vec3 V)
{
	Mat3 Vh;
	Vh<<0, -V(2), V(1),
	   V(2),   0,-V(0),
	   -V(1), V(0),  0;
	return Vh;
}

Vec3 SO3::vee(Mat3 Vh)
{
	Vec3 V;
	V<<Vh(2,1),Vh(0,2),Vh(1,0);
	return V;
}

Mat4 SO3::toSE3(Mat3 g)
{
	Mat4 gSE3=Mat4::Zero();
	gSE3(3,3)=1;
	gSE3.topLeftCorner(3,3)=g;
	return gSE3;
}

Vec6 SO3::tose3(Vec3 V)
{
	Vec6 v;
	v<<0,0,0,V;
	return v;
}

Mat3 SO3::inverse(Mat3 g)
{
	return g.transpose();
}

Vec3 SO3::getRPY(Mat3 g)
{
	Vec3 w=Vec3::Zero();
	double TR=g.trace();
	if(fabs(TR-3)>1e-10) 
	{
		double th=acos((TR-1)/2);
		w(0)=g(2,1)-g(1,2);
		w(1)=g(0,2)-g(2,0);
		w(2)=g(1,0)-g(0,1);
		w=w*th/(2*sin(th));
	}
	return w;
}

class SE3
{
	public:
		static const uchar row;
		static const uchar col;
		static const uchar dim;

		static const Vec6 e[6]; 
		static const Mat4 E[6]; 


		static Mat4 exp(Vec6 v);
		static Vec6 log(Mat4 g);
		static Mat6 dexp(Vec6 v);
		static Mat6 dexpinv(Vec6 v);
		static Mat6 dexpinvAd(Vec6 v);
		static Eigen::Matrix<double,6,36> ddexp(Vec6 xi);
		static Eigen::Matrix<double,6,36> ddexpinv(Vec6 xi);

		static Mat6 ad(Vec6 v);
		static Mat6 Ad(Mat4 g);

		static Mat4 hat(Vec6 v);
		static Vec6 vee(Mat4 Vh);

		static Mat4 cay(Vec6 v);
		static Vec6 cayinv(Mat4 g);
		static Mat6 dcay(Vec6 v);
		static Mat6 dcayinv(Vec6 v);

		static Mat4 inverse(Mat4 g);

};

const uchar SE3::row=4;
const uchar SE3::col=4;
const uchar SE3::dim=6;

const Vec6 SE3::e[6]={(Vec6()<<1,0,0,0,0,0).finished(),
				 (Vec6()<<0,1,0,0,0,0).finished(),
				 (Vec6()<<0,0,1,0,0,0).finished(),
				 (Vec6()<<0,0,0,1,0,0).finished(),
				 (Vec6()<<0,0,0,0,1,0).finished(),
				 (Vec6()<<0,0,0,0,0,1).finished()
				};

const Mat4 SE3::E[6]={(Mat4()<< 0,  0, 0, 0,
							    0,  0,-1, 0,
							    0,  1, 0, 0,
								0,  0, 0, 0).finished(),

					  (Mat4()<< 0,  0, 1, 0,
							    0,  0, 0, 0,
							   -1,  0, 0, 0,
								0,  0, 0, 0).finished(),

					  (Mat4()<< 0, -1, 0, 0,
							    1,  0, 0, 0,
								0,  0, 0, 0,
								0,  0, 0, 0).finished(),

					  (Mat4()<< 0,  0, 0, 1,
							    0,  0, 0, 0,
								0,  0, 0, 0,
								0,  0, 0, 0).finished(),

					  (Mat4()<< 0,  0, 0, 0,
							    0,  0, 0, 1,
								0,  0, 0, 0,
								0,  0, 0, 0).finished(),

					  (Mat4()<< 0,  0, 0, 0,
							    0,  0, 0, 0,
								0,  0, 0, 1,
								0,  0, 0, 0).finished()
                };

Mat4 SE3::exp(Vec6 v)
{
	Mat4 g=Mat4::Identity();

	double th=v.head(3).norm();

	if(th<1e-10)
		g.block(0,3,3,1)=v.tail(3);
	else
	{
		Vec3 w=v.head(3);
		g.block(0,0,3,3)=SO3::exp(w);
		g.block(0,3,3,1)=(((Mat3::Identity()-g.block(0,0,3,3))*SO3::hat(w)+ w*w.transpose())/(th*th))*v.tail(3);
	}
	return g;
}

Vec6 SE3::log(Mat4 g)
{
	Vec6 v;
	v.head(3)=SO3::log(g.block(0,0,3,3));

	double th=v.head(3).norm();
	if(th<1e-7)
		v.tail(3)=g.block(0,3,3,1);
	else
	{
		Vec3 w=v.head(3);
		v.tail(3)=(((Mat3::Identity()-g.block(0,0,3,3))*SO3::hat(w)+w*w.transpose())/(th*th)).inverse()*g.block(0,3,3,1);
	}

	return v;
}

Mat6 SE3::dexp(Vec6 xi)
{
	Vec3 w=xi.head(3);
	Mat3 W=SO3::hat(w);
	Mat3 V=SO3::hat(xi.tail(3));

	double th=w.norm();

	if(th<1e-6)
		return Mat6::Identity()+ad(xi)/2;
	else
	{
		Mat3 W2=W*W;

		Mat3 WV=W*V;
		Mat3 VW=V*W;
		Mat3 W2V=W2*V;
		Mat3 WVW=W*VW;
		Mat3 VW2=V*W2;
		Mat3 WVW2=W*VW2;
		Mat3 W2VW=W2*VW;

		double sth=sin(th);
		double cth=cos(th);

		double th2=th*th;
		double th3=th2*th;
		double th4=th3*th;
		double th5=th4*th;
		double th6=th5*th;

		Mat6 DEXP=Mat6::Zero();
		
		DEXP.block(0,0,3,3)=Mat3::Identity()+((1-cth)/th2)*W+((th-sth)/th3)*W2;
		DEXP.block(3,3,3,3)=DEXP.block(0,0,3,3);

		double c1=(2.0-2.0*cth-th*sth/2.0)/th2;
		double c2=(th-sth)/th3;
		double c3=(1-cth-th*sth/2.0)/th4;
		double c4=(th-3.0*sth/2.0+th*cth/2.0)/th5;

		DEXP.block(3,0,3,3)=c1*V+
							c2*(WV+VW)+
							c3*(W2V+WVW+VW2)+
							c4*(W2VW+WVW2);

		return DEXP;
	}
}

Mat6 SE3::dexpinv(Vec6 xi)
{
	Vec3 w=xi.head(3);
	Mat3 W=SO3::hat(w);
	Mat3 V=SO3::hat(xi.tail(3));

	double th=w.norm();

	if(th<1e-6)
		return Mat6::Identity()-ad(xi)/2;
	else
	{
		Mat3 W2=W*W;

		Mat3 WV=W*V;
		Mat3 VW=V*W;
		Mat3 W2VW=W2*VW;
		Mat3 WVW2=WV*W2;


		double sth=sin(th);
		double cth=cos(th);
		
		double th2=th*th;
		double th3=th2*th;
		double th4=th3*th;

		Mat6 DEXPINV=Mat6::Zero();

		DEXPINV.block(0,0,3,3)=Mat3::Identity()-W/2+((th*sth+2*cos(th)-2)/(2*th2*(cos(th)-1)))*W2;
		DEXPINV.block(3,3,3,3)=DEXPINV.block(0,0,3,3);

		DEXPINV.block(3,0,3,3)=-V/2
							   +((th*sth+2*cth-2)/(2*th2*(cth-1)))*(WV+VW)
							   +((th2+th*sth+4*cth-4)/(4*th4*(cth-1)))*(W2VW+WVW2);

		return DEXPINV;
	}
}

Mat6 SE3::dexpinvAd(Vec6 v)
{
	return dexpinv(-v);
}

Eigen::Matrix<double,6,36> SE3::ddexp(Vec6 xi)
{
	Vec3 w=xi.head(3);
	Vec3 v=xi.tail(3);

	double th=w.norm();
	double th2=th*th;
	double th3=th2*th;
	double th4=th3*th;
	double th5=th4*th;
	double th7=th5*th2;

	Mat3 W=SO3::hat(w);
	Mat3 V=SO3::hat(v);

	double c1=0;
	double c2=0;
	double c3=0;
	double c4=0;

	double dc1=0;
	double dc2=0;
	double dc3=0;
	double dc4=0;

	double sth=sin(th);
	double cth=cos(th);

	Eigen::Matrix<double,6,36> DF=Eigen::Matrix<double,6,36>::Zero();

	Eigen::Matrix<double,3,9> ddexp_SO3=SO3::ddexp(w);

	DF.block(0,0,3,3)=DF.block(3,3,3,3)=ddexp_SO3.block(0,0,3,3);
	DF.block(0,6,3,3)=DF.block(3,9,3,3)=ddexp_SO3.block(0,3,3,3);
	DF.block(0,12,3,3)=DF.block(3,15,3,3)=ddexp_SO3.block(0,6,3,3);

	if (th>1e-6)
	{
		c1=(2.0-2.0*cth-th*sth/2.0)/th2;
		c2=(th-sth)/th3;
		c3=(1-cth-th*sth/2.0)/th4;
		c4=(th-3.0*sth/2.0+th*cth/2.0)/th5;

		dc1=(8*cth-th2*cth+5*th*sth-8)/(2*th4);
		dc2=-(th*(cth-1)+3*(th-sth))/th5;
		dc3=dc1/th2;
		dc4=-(8*th-15*sth+th2*sth+7*th*cth)/(2*th7);
	
		Mat3 W2=W*W;

		Mat3 WV=W*V;
		Mat3 VW=V*W;
		Mat3 W2V=W2*V;
		Mat3 WVW=W*VW;
		Mat3 VW2=V*W2;
		Mat3 WVW2=W*VW2;
		Mat3 W2VW=W2*VW;


		Mat3 DTH=dc1*V+dc2*(WV+VW)+dc3*(W2V+WVW+VW2)+dc4*(W2VW+WVW2);
		
		Mat3 DF1=c2*(SO3::E[0]*V+V*SO3::E[0])+
				 c3*(SO3::E[0]*(WV+VW)+(WV+VW)*SO3::E[0]+W*SO3::E[0]*V+V*SO3::E[0]*W)+
				 c4*(SO3::E[0]*(WVW+VW2)+(W2V+WVW)*SO3::E[0]+W*(SO3::E[0]*V+V*SO3::E[0])*W)+
				 DTH*w(0);

		Mat3 DF2=c2*(SO3::E[1]*V+V*SO3::E[1])+
				 c3*(SO3::E[1]*(WV+VW)+(WV+VW)*SO3::E[1]+W*SO3::E[1]*V+V*SO3::E[1]*W)+
				 c4*(SO3::E[1]*(WVW+VW2)+(W2V+WVW)*SO3::E[1]+W*(SO3::E[1]*V+V*SO3::E[1])*W)+
				 DTH*w(1);

		Mat3 DF3=c2*(SO3::E[2]*V+V*SO3::E[2])+
				 c3*(SO3::E[2]*(WV+VW)+(WV+VW)*SO3::E[2]+W*SO3::E[2]*V+V*SO3::E[2]*W)+
				 c4*(SO3::E[2]*(WVW+VW2)+(W2V+WVW)*SO3::E[2]+W*(SO3::E[2]*V+V*SO3::E[2])*W)+
				 DTH*w(2);

		Mat3 DF4=c1*SO3::E[0]+
				 c2*(W*SO3::E[0]+SO3::E[0]*W)+
				 c3*(W2*SO3::E[0]+W*SO3::E[0]*W+SO3::E[0]*W2)+
				 c4*(W2*SO3::E[0]*W+W*SO3::E[0]*W2);

		Mat3 DF5=c1*SO3::E[1]+
				 c2*(W*SO3::E[1]+SO3::E[1]*W)+
				 c3*(W2*SO3::E[1]+W*SO3::E[1]*W+SO3::E[1]*W2)+
				 c4*(W2*SO3::E[1]*W+W*SO3::E[1]*W2);

		Mat3 DF6=c1*SO3::E[2]+
				 c2*(W*SO3::E[2]+SO3::E[2]*W)+
				 c3*(W2*SO3::E[2]+W*SO3::E[2]*W+SO3::E[2]*W2)+
				 c4*(W2*SO3::E[2]*W+W*SO3::E[2]*W2);

		DF.block(3,0,3,3)=DF1;
		DF.block(3,6,3,3)=DF2;
		DF.block(3,12,3,3)=DF3;
		DF.block(3,18,3,3)=DF4;
		DF.block(3,24,3,3)=DF5;
		DF.block(3,30,3,3)=DF6;
	}
	else
	{
		c1=1.0/2.0;
		c2=1.0/6.0;

		Mat3 DF1=c2*(SO3::E[0]*V+V*SO3::E[0]);
		Mat3 DF2=c2*(SO3::E[1]*V+V*SO3::E[1]);
		Mat3 DF3=c2*(SO3::E[2]*V+V*SO3::E[2]);
		Mat3 DF4=c1*SO3::E[0];
		Mat3 DF5=c1*SO3::E[1];
		Mat3 DF6=c1*SO3::E[2];

		DF.block(3,0,3,3)=DF1;
		DF.block(3,6,3,3)=DF2;
		DF.block(3,12,3,3)=DF3;
		DF.block(3,18,3,3)=DF4;
		DF.block(3,24,3,3)=DF5;
		DF.block(3,30,3,3)=DF6;
	}

	return DF;
}

Eigen::Matrix<double,6,36> SE3::ddexpinv(Vec6 xi)
{
	Vec3 w=xi.head(3);
	Vec3 v=xi.tail(3);

	double th=w.norm();
	double th2=th*th;
	double th3=th2*th;
	double th4=th3*th;
	double th6=th4*th2;

	Mat3 W=SO3::hat(w);
	Mat3 V=SO3::hat(v);

	double c1=0;
	double c2=0;
	double c3=0;

	double dc1=0;
	double dc2=0;

	Eigen::Matrix<double,6,36> DF=Eigen::Matrix<double,6,36>::Zero();

	Eigen::Matrix<double,3,9> ddexpinv_SO3=SO3::ddexpinv(w);

	DF.block(0,0,3,3)=DF.block(3,3,3,3)=ddexpinv_SO3.block(0,0,3,3);
	DF.block(0,6,3,3)=DF.block(3,9,3,3)=ddexpinv_SO3.block(0,3,3,3);
	DF.block(0,12,3,3)=DF.block(3,15,3,3)=ddexpinv_SO3.block(0,6,3,3);

	if (th>1e-6)
	{
		double sth=sin(th);
		double cth=cos(th)-1;

		c1=-0.5;
		c2=(th*sth+2*cth)/(2*th2*cth);
		c3=(th2+th*sth+4*cth)/(4*th4*cth);

		dc1=-(4*cth+th*sth+th2)/(2*th4*cth);
		dc2=(th3*sth-3*th*(th+sth)*cth-16*cth*cth)/(4*th6*cth*cth);

		Mat3 WV=W*V;
		Mat3 VW=V*W;
		Mat3 W2=W*W;
		Mat3 W2V=W*WV;
		Mat3 WVW=WV*W;
		Mat3 VW2=VW*W;
		Mat3 W2VW=W2V*W;
		Mat3 WVW2=W*VW2;
		Mat3 DTH=dc1*(WV+VW)+dc2*(W2VW+WVW2);

		Mat3 DF1=c2*(SO3::E[0]*V+V*SO3::E[0])+
				 c3*(SO3::E[0]*(WVW+VW2)+(W2V+WVW)*SO3::E[0]+W*(SO3::E[0]*V+V*SO3::E[0])*W)+
				 DTH*w(0);

		Mat3 DF2=c2*(SO3::E[1]*V+V*SO3::E[1])+
				 c3*(SO3::E[1]*(WVW+VW2)+(W2V+WVW)*SO3::E[1]+W*(SO3::E[1]*V+V*SO3::E[1])*W)+
				 DTH*w(1);

		Mat3 DF3=c2*(SO3::E[2]*V+V*SO3::E[2])+
				 c3*(SO3::E[2]*(WVW+VW2)+(W2V+WVW)*SO3::E[2]+W*(SO3::E[2]*V+V*SO3::E[2])*W)+
				 DTH*w(2);

		Mat3 DF4=c1*SO3::E[0]+
				 c2*(W*SO3::E[0]+SO3::E[0]*W)+
				 c3*(W2*SO3::E[0]*W+W*SO3::E[0]*W2);

		Mat3 DF5=c1*SO3::E[1]+
				 c2*(W*SO3::E[1]+SO3::E[1]*W)+
				 c3*(W2*SO3::E[1]*W+W*SO3::E[1]*W2);

		Mat3 DF6=c1*SO3::E[2]+
				 c2*(W*SO3::E[2]+SO3::E[2]*W)+
				 c3*(W2*SO3::E[2]*W+W*SO3::E[2]*W2);

		DF.block(3,0,3,3)=DF1;
		DF.block(3,6,3,3)=DF2;
		DF.block(3,12,3,3)=DF3;
		DF.block(3,18,3,3)=DF4;
		DF.block(3,24,3,3)=DF5;
		DF.block(3,30,3,3)=DF6;
	}
	else
	{
		c1=-1.0/2.0;
		c2=1.0/12.0;

		Mat3 DF1=c2*(SO3::E[0]*V+V*SO3::E[0]);
		Mat3 DF2=c2*(SO3::E[1]*V+V*SO3::E[1]);
		Mat3 DF3=c2*(SO3::E[2]*V+V*SO3::E[2]);
		Mat3 DF4=c1*SO3::E[0];
		Mat3 DF5=c1*SO3::E[1];
		Mat3 DF6=c1*SO3::E[2];

		DF.block(3,0,3,3)=DF1;
		DF.block(3,6,3,3)=DF2;
		DF.block(3,12,3,3)=DF3;
		DF.block(3,18,3,3)=DF4;
		DF.block(3,24,3,3)=DF5;
		DF.block(3,30,3,3)=DF6;
	}

	return DF;
}


Mat6 SE3::ad(Vec6 v)
{
	return (Mat6()<<SO3::hat(v.head(3)),
					Mat3::Zero(),
					SO3::hat(v.tail(3)),
					SO3::hat(v.head(3))).finished();
}

Mat6 SE3::Ad(Mat4 g)
{
	Mat3 R=g.block(0,0,3,3);
	Vec3 p=g.block(0,3,3,1);

	return (Mat6()<<R,
					Mat3::Zero(),
					SO3::hat(p)*R,
					R).finished();
}

Mat4 SE3::hat(Vec6 v)
{
	return (Mat4()<<SO3::hat(v.head(3)), v.tail(3),
					0, 0, 0, 0 ).finished();
}

Vec6 SE3::vee(Mat4 Vh)
{
	return (Vec6()<<SO3::vee(Vh.block(0,0,3,3)), Vh.block(0,3,3,1)).finished();
}

Mat4 SE3::cay(Vec6 v)
{
	Vec3 w=v.head(3);
	double th=w.squaredNorm();

	return (Mat4()<< SO3::cay(w),
					 (Mat3::Identity()+SO3::ad(w)/2+w*w.transpose()/4)*v.tail(3)*4/(4+th),
					 0, 0, 0, 1).finished();
}

Vec6 SE3::cayinv(Mat4 g)
{
	return SE3::vee(-2*(Mat4::Identity()+g).inverse()*(Mat4::Identity()-g)); 
}

Mat6 SE3::dcay(Vec6 v)
{
	double th=(v.head(3).squaredNorm()+4);

	Vec3 w=v.head(3);
	Mat3 vh=SO3::hat(v.tail(3));
	Mat3 A=Mat3::Identity()+SO3::hat(v.head(3))/2;

	Mat6 F=Mat6::Zero();

	F.block(0,0,3,3)=A*(4/th); 
	F.block(3,0,3,3)=vh*A*(2/th);
	F.block(3,3,3,3)=(w*w.transpose()+4*A)/th;

	return F;
}

Mat6 SE3::dcayinv(Vec6 v)
{
	Vec3 w=v.head(3);
	Mat3 vh=SO3::hat(v.tail(3));
	Mat3 A=Mat3::Identity()-SO3::hat(v.head(3))/2;

	Mat6 F=Mat6::Zero();
	F.block(0,0,3,3)=A+w*w.transpose()/4;
	F.block(3,0,3,3)=-A*vh/2;
	F.block(3,3,3,3)=A;

	return F;
}

Mat4 SE3::inverse(Mat4 g)
{
	Mat3 RT=g.block(0,0,3,3).transpose();
	return (Mat4()<<RT, -RT*g.block(0,3,3,1), 0, 0, 0, 1).finished();
}
#endif
