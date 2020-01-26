#ifndef ASFEM_MATHUTILS_H
#define ASFEM_MATHUTILS_H


#include <iostream>
#include <iomanip>
#include <cmath>

#include "RankFourTensor.h"

#include "Eigen/Eigen"

using namespace std;

inline double operator*(const Eigen::Vector3d &a,const Eigen::Vector3d &b){
    return a.coeff(0)*b.coeff(0)+a.coeff(1)*b.coeff(1)+a.coeff(2)*b.coeff(2);
}

inline Eigen::Vector3d operator*(const Eigen::Vector3d &a,const double &b){
    return a*b;
}

inline Eigen::Vector3d operator*(const double &a,const Eigen::Vector3d &b){
    return b*a;
}

inline Eigen::VectorXd operator*(const Eigen::VectorXd &a,const Eigen::VectorXd &b){
    Eigen::VectorXd temp=a;
    for(unsigned int i=0;i<a.size();++i){
        temp.coeffRef(i)*=b.coeff(i);
    }
    return temp;
}

inline double AtDotB(const Eigen::VectorXd &a,const Eigen::VectorXd &b){
    double temp=0.0;
    for(unsigned int i=0;i<a.size();++i){
        temp+=a.coeff(i)*b.coeff(i);
    }
    return temp;
}

inline double GetIKjlComponent(const RankFourTensor &jac,
                               const int &i,const int &k,
                               const Eigen::Vector3d &grad_test,
                               const Eigen::Vector3d &grad_phi){
    return (jac(i,1,k,1)*grad_phi.coeff(0)
           +jac(i,1,k,2)*grad_phi.coeff(1)
           +jac(i,1,k,3)*grad_phi.coeff(2))*grad_test.coeff(0)
          +(jac(i,2,k,1)*grad_phi.coeff(0)
           +jac(i,2,k,2)*grad_phi.coeff(1)
           +jac(i,2,k,3)*grad_phi.coeff(2))*grad_test.coeff(1)
          +(jac(i,3,k,1)*grad_phi.coeff(0)
           +jac(i,3,k,2)*grad_phi.coeff(1)
           +jac(i,3,k,3)*grad_phi.coeff(2))*grad_test.coeff(2);
}

#endif // ASFEM_MATHUTILS_H