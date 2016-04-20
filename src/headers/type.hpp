/*************************************************************************
    > File Name: type.hpp
    > Author: Taosha Fan
    > Mail: tfan1@jhu.edu 
    > Created Time: Fri 20 Feb 2015 08:54:49 PM EST
 ************************************************************************/
#ifndef TYPE
#define TYPE

#include<Eigen/Dense>

typedef Eigen::MatrixXd Mat;
typedef Eigen::Matrix<double,1,1> Mat1;
typedef Eigen::Matrix2d Mat2;
typedef Eigen::Matrix3d Mat3;
typedef Eigen::Matrix4d Mat4;
typedef Eigen::Matrix<double,6,6> Mat6;
typedef Eigen::Matrix<double,9,9> Mat9;
typedef Eigen::Matrix<double,12,12> Mat12;

typedef Eigen::VectorXd Vec;
typedef Eigen::Matrix<double,1,1> Vec1;
typedef Eigen::Vector2d Vec2;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector4d Vec4;
typedef Eigen::Matrix<double,6,1> Vec6;
typedef Eigen::Matrix<double,9,1> Vec9;
typedef Eigen::Matrix<double,12,1> Vec12;

typedef unsigned char uchar;
#endif
