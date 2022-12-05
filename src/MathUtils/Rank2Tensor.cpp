//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.10.17
//+++ Update : 2022.07.24->re-write & re-design @Yang Bai
//+++ Purpose: Implement rank-2 tensor class for the common
//+++          tensor manipulation in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MathUtils/Rank2Tensor.h"
#include "Eigen/Eigen"

Rank2Tensor::Rank2Tensor(){
    m_vals.resize(9,0.0);
}
Rank2Tensor::Rank2Tensor(const double &val){
    m_vals.resize(9,val);
}
Rank2Tensor::Rank2Tensor(const Rank2Tensor &a){
    m_vals.resize(9,0.0);
    for(int i=0;i<N2;i++) m_vals[i]=a.m_vals[i];
}
Rank2Tensor::Rank2Tensor(const InitMethod &initmethod){
    if(initmethod==InitMethod::ZERO){
        m_vals.resize(9,0.0);
    }
    else if(initmethod==InitMethod::IDENTITY){
        m_vals.resize(9,0.0);
        setToIdentity();
    }
    else if(initmethod==InitMethod::RANDOM){
        m_vals.resize(9,0.0);
        setToRandom();
    }
    else{
        MessagePrinter::printErrorTxt("unsupported initialize method in rank-2 tensor");
        MessagePrinter::exitAsFem();
    }
}
Rank2Tensor::~Rank2Tensor(){
    m_vals.clear();
}
//**********************************************************************
Rank2Tensor operator*(const double &lhs,const Rank2Tensor &a){
    Rank2Tensor temp(0.0);
    for(int i=0;i<a.N2;i++) temp.m_vals[i]=lhs*a.m_vals[i];
    return temp;
}
Vector3d operator*(const Vector3d &lhs,const Rank2Tensor &a){
    Vector3d temp(0.0);
    for(int j=1;j<=a.N;j++){
        temp(j)=lhs(1)*a(1,j)+lhs(2)*a(2,j)+lhs(3)*a(3,j);
    }
    return temp;
}
Rank2Tensor Rank2Tensor::doubledot(const Rank4Tensor &a) const{
        // return A:B calculation
        Rank2Tensor temp(0.0);
        for(int i=1;i<=N;i++){
            for(int j=1;j<=N;j++){
                for(int k=1;k<=N;k++){
                    for(int l=1;l<=N;l++){
                        temp(k,l)+=(*this)(i,j)*a(i,j,k,l);
                    }
                }
            }
        }
        return temp;
    }
//************************************************************************
void Rank2Tensor::setFromGradU(const Vector3d &gradUx){
    (*this)(1,1)=gradUx(1);(*this)(1,2)=0.0;(*this)(1,3)=0.0;
    (*this)(2,1)=      0.0;(*this)(2,2)=0.0;(*this)(2,3)=0.0;
    (*this)(3,1)=      0.0;(*this)(3,2)=0.0;(*this)(3,3)=0.0;
}
void Rank2Tensor::setFromGradU(const Vector3d &gradUx,const Vector3d &gradUy){
    (*this)(1,1)=gradUx(1);(*this)(1,2)=gradUx(2);(*this)(1,3)=0.0;
    (*this)(2,1)=gradUy(1);(*this)(2,2)=gradUy(2);(*this)(2,3)=0.0;
    (*this)(3,1)=      0.0;(*this)(3,2)=      0.0;(*this)(3,3)=0.0;
}
void Rank2Tensor::setFromGradU(const Vector3d &gradUx,const Vector3d &gradUy,const Vector3d &gradUz){
    (*this)(1,1)=gradUx(1);(*this)(1,2)=gradUx(2);(*this)(1,3)=gradUx(3);
    (*this)(2,1)=gradUy(1);(*this)(2,2)=gradUy(2);(*this)(2,3)=gradUy(3);
    (*this)(3,1)=gradUz(1);(*this)(3,2)=gradUz(2);(*this)(3,3)=gradUz(3);
}
void Rank2Tensor::setRotationTensorFromEulerAngle(const double &theta1,const double &theta2,const double &theta3){
    // set a rotation tensor from Euler-angle
    const double PI=3.14159265359;
    double x1=cos(theta1*PI/180.0);
    double x2=cos(theta2*PI/180.0);
    double x3=cos(theta3*PI/180.0);
    double y1=sin(theta1*PI/180.0);
    double y2=sin(theta2*PI/180.0);
    double y3=sin(theta3*PI/180.0);

    (*this)(1,1)= x1*x3-x2*y1*y3;
    (*this)(1,2)= x3*y1+x1*x2*y3;
    (*this)(1,3)= y2*y3;

    (*this)(2,1)=-x1*y3-x2*x3*y1;
    (*this)(2,2)= x1*x2*x3-y1*y3;
    (*this)(2,3)= x3*y2;

    (*this)(3,1)= y1*y2;
    (*this)(3,2)=-x1*y2;
    (*this)(3,3)= x2;
}
//*******************************************************************
//*** for advanced math operators
//*******************************************************************
Rank2Tensor exp(const Rank2Tensor &a){
    Rank2Tensor I;
    I.setToIdentity();
    return I
          +a
          +a*a*(1.0/(1.0*2.0))
          +a*a*a*(1.0/(1.0*2.0*3.0))
          +a*a*a*a*(1.0/(1.0*2.0*3.0*4.0))
          +a*a*a*a*a*(1.0/(1.0*2.0*3.0*4.0*5.0))
          +a*a*a*a*a*a*(1.0/(1.0*2.0*3.0*4.0*5.0*6.0));
}
Rank2Tensor dexp(const double &a,const Rank2Tensor &b){
    // return dexp(ab)/db
    return b
          +b*b*a*(1.0/1.0)
          +b*b*b*a*a*(1.0/(1.0*2.0))
          +b*b*b*b*a*a*a*(1.0/(1.0*2.0*3.0))
          +b*b*b*b*b*a*a*a*a*(1.0/(1.0*2.0*3.0*4.0))
          +b*b*b*b+b*b*a*a*a*a*a*(1.0/(1.0*2.0*3.0*4.0*5.0));
}
//*******************************************************************
//*** some higher order tensor calculations
//*******************************************************************
Rank4Tensor Rank2Tensor::otimes(const Rank2Tensor &a) const{
    // return C_ijkl=a_ij*b_kl
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(i,j)*a(k,l);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank2Tensor::ijXlk(const Rank2Tensor &a) const{
    // return C_ijkl=a_ij*b_lk
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(i,j)*a(l,k);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank2Tensor::odot(const Rank2Tensor &a) const{
    // extremely useful for:
    //    dA^-1/dA=-A^-1 \otimes A^-1
    //            =-0.5*Ainv\otimes Ainv=rank-4 tensor
    // for nonlinear constitutive law
    // the proof can be found here:
    // https://en.wikipedia.org/wiki/Tensor_derivative_(continuum_mechanics)
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=0.5*((*this)(i,k)*a(j,l)+(*this)(i,l)*a(j,k));
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For IkJl and IklJ
//********************************************************
Rank4Tensor Rank2Tensor::ikXjl(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(i,k)*a(j,l);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank2Tensor::ikXlj(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(i,k)*a(l,j);
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For IlJk and IlkJ
//********************************************************
Rank4Tensor Rank2Tensor::ilXjk(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(i,l)*a(j,k);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank2Tensor::ilXkj(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(i,l)*a(k,j);
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For JIkl and JIlk
//********************************************************
Rank4Tensor Rank2Tensor::jiXkl(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(j,i)*a(k,l);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank2Tensor::jiXlk(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(j,i)*a(l,k);
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For JkIl and JklI
//********************************************************
Rank4Tensor Rank2Tensor::jkXil(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(j,k)*a(i,l);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank2Tensor::jkXli(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(j,k)*a(l,i);
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For JlIk and JlkI
//********************************************************
Rank4Tensor Rank2Tensor::jlXik(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(j,l)*a(i,k);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank2Tensor::jlXki(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(j,l)*a(k,i);
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For kIJl and kIlJ
//********************************************************
Rank4Tensor Rank2Tensor::kiXjl(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(k,i)*a(j,l);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank2Tensor::kiXlj(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(k,i)*a(l,j);
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For kJIl and kJlI
//********************************************************
Rank4Tensor Rank2Tensor::kjXil(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(k,j)*a(i,l);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank2Tensor::kjXli(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(k,j)*a(l,i);
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For klIJ and klJI
//********************************************************
Rank4Tensor Rank2Tensor::klXij(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(k,l)*a(i,j);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank2Tensor::klXji(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(k,l)*a(j,i);
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For lIkJ and lIJk
//********************************************************
Rank4Tensor Rank2Tensor::liXkj(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(l,i)*a(k,j);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank2Tensor::liXjk(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(l,i)*a(j,k);
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For lJIk and lJkI
//********************************************************
Rank4Tensor Rank2Tensor::ljXik(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(l,j)*a(i,k);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank2Tensor::ljXki(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(l,j)*a(k,i);
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For lkIJ and lkJI
//********************************************************
Rank4Tensor Rank2Tensor::lkXij(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(l,k)*a(i,j);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank2Tensor::lkXji(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=(*this)(l,k)*a(j,i);
                }
            }
        }
    }
    return temp;
}
//**************************************************************
//*** For eigen value and eigen vectors and other 
//*** stress and strain decomposition related functions
//**************************************************************
void Rank2Tensor::calcEigenValueAndEigenVectors(double (&eigval)[3],Rank2Tensor &eigvec) const{
    Eigen::Matrix3d _M;

    _M<<(*this)(1,1),(*this)(1,2),(*this)(1,3),
        (*this)(2,1),(*this)(2,2),(*this)(2,3),
        (*this)(3,1),(*this)(3,2),(*this)(3,3);
    
    Eigen::EigenSolver<Eigen::Matrix3d> _eigen_solver;
    _eigen_solver.compute(_M);
    //
    eigval[0]=_eigen_solver.eigenvalues()(0).real();
    eigval[1]=_eigen_solver.eigenvalues()(1).real();
    eigval[2]=_eigen_solver.eigenvalues()(2).real();
    for(int i=0;i<N;i++){
        eigvec(1,i+1)=_eigen_solver.eigenvectors()(0,i).real();
        eigvec(2,i+1)=_eigen_solver.eigenvectors()(1,i).real();
        eigvec(3,i+1)=_eigen_solver.eigenvectors()(2,i).real();
    }
}
Rank4Tensor Rank2Tensor::calcPositiveProjTensor(double (&eigval)[3],Rank2Tensor &eigvec) const{
    // Algorithm is taken from:
    // C. Miehe and M. Lambrecht, Commun. Numer. Meth. Engng 2001; 17:337~353
    // https://onlinelibrary.wiley.com/doi/epdf/10.1002/cnm.404

    calcEigenValueAndEigenVectors(eigval,eigvec);

    // C=F^T F=lambda_i M_i where M_i=n_i x n_i
    //         lambda_i-->eigenvalue  n_i--> the eigenvector

    double epos[3],diag[3];
    for(int i=0;i<N;i++){
        epos[i]=0.5*(abs(eigval[i])+eigval[i]);
        diag[i]=0.0;
        if(eigval[i]>0.0){
            diag[i]=1.0;
        }
    }
    Rank4Tensor ProjPos(0.0);
    Rank2Tensor Ma(0.0),Mb(0.0);

    ProjPos.setToZeros();

    // calculate Ma defined in Eq.(9)-2
    // Ma=n_a x n_a
    for(int i=1;i<=N;i++){
        Ma.setFromVectorDyad(eigvec.getIthCol(i),eigvec.getIthCol(i));
        // Eq.(19), first term on the right side
        ProjPos+=Ma.otimes(Ma)*diag[i-1];
    }

    // Now we calculate the Gab and Gba
    // We need a new rank-2 tensor Mb(same defination as Ma)
    double theta_ab;// defined in Eq.(21)-1
    Rank4Tensor Gab(0.0),Gba(0.0);
    const double tol=1.0e-13;
    for(int a=0;a<N;a++){
        for(int b=0;b<a;b++){
            Ma.setFromVectorDyad(eigvec.getIthCol(a+1),eigvec.getIthCol(a+1));// Eq.(12)
            Mb.setFromVectorDyad(eigvec.getIthCol(b+1),eigvec.getIthCol(b+1));// change the order of Eq.(12)

            Gab=Ma.ikXjl(Mb)+Ma.ilXjk(Mb);
            Gba=Mb.ikXjl(Ma)+Mb.ilXjk(Ma);
            // since only positive term is involved
            // e_a=0.5*(abs(lambda_a)+lambda_a)
            // P is defined as: 2dE/dC in Eq.(8)-2
            //  but 2dM/dC=(Gab+Gba)/(lambda_a-lambda_b)
            // E(C)=sum(e_a*M_a)
            // P=2dE(C)/dC=2(dE(C)/dM)*(dM/dC)
            if(abs(eigval[a]-eigval[b])<=tol){
                //if limit lambda_a to lambda_b in Eq.(24)
                theta_ab=0.5*(diag[a]+diag[b])/2.0;
            }
            else{
                // ea(lambda_a)=(1/2)*(lambda_a-1)
                // m=2 for green strain case in Eq.(16)
                theta_ab=0.5*(epos[a]-epos[b])/(eigval[a]-eigval[b]);// Eq.(21)-1
            }
            ProjPos+=theta_ab*(Gab+Gba);
        }
    }
    return ProjPos;
}

Rank4Tensor Rank2Tensor::getPositiveProjectionTensor() const{
    // Algorithm is taken from:
    // C. Miehe and M. Lambrecht, Commun. Numer. Meth. Engng 2001; 17:337~353
    // https://onlinelibrary.wiley.com/doi/epdf/10.1002/cnm.404

    double eigval[3];
    Rank2Tensor eigvec;

    calcEigenValueAndEigenVectors(eigval,eigvec);
    // C=F^T F=lambda_i M_i where M_i=n_i x n_i
    //         lambda_i-->eigenvalue  
    //         n_i     --> the eigenvector
    double epos[3],diag[3];
    for(int i=0;i<N;i++){
        epos[i]=0.5*(abs(eigval[i])+eigval[i]);
        diag[i]=0.0;
        if(eigval[i]>0.0){
            diag[i]=1.0;
        }
    }
    Rank4Tensor ProjPos(0.0);
    Rank2Tensor Ma(0.0),Mb(0.0);

    ProjPos.setToZeros();

    // calculate Ma defined in Eq.(9)-2
    // Ma=n_a x n_a
    for(int i=1;i<=N;i++){
        Ma.setFromVectorDyad(eigvec.getIthCol(i),eigvec.getIthCol(i));
        // Eq.(19), first term on the right side
        ProjPos+=Ma.otimes(Ma)*diag[i-1];
    }

    // Now we calculate the Gab and Gba
    // We need a new rank-2 tensor Mb(same defination as Ma)
    double theta_ab;// defined in Eq.(21)-1
    Rank4Tensor Gab(0.0),Gba(0.0);
    const double tol=1.0e-13;
    for(int a=0;a<N;a++){
        for(int b=0;b<a;b++){
            Ma.setFromVectorDyad(eigvec.getIthCol(a+1),eigvec.getIthCol(a+1));// Eq.(12)
            Mb.setFromVectorDyad(eigvec.getIthCol(b+1),eigvec.getIthCol(b+1));// change the order of Eq.(12)

            Gab=Ma.ikXjl(Mb)+Ma.ilXjk(Mb);
            Gba=Mb.ikXjl(Ma)+Mb.ilXjk(Ma);
            // since only positive term is involved
            // e_a=0.5*(abs(lambda_a)+lambda_a)
            // P is defined as: 2dE/dC in Eq.(8)-2
            //  but 2dM/dC=(Gab+Gba)/(lambda_a-lambda_b)
            // E(C)=sum(e_a*M_a)
            // P=2dE(C)/dC=2(dE(C)/dM)*(dM/dC)
            if(abs(eigval[a]-eigval[b])<=tol){
                //if limit lambda_a to lambda_b in Eq.(24)
                theta_ab=0.5*(diag[a]+diag[b])/2.0;
            }
            else{
                // ea(lambda_a)=(1/2)*(lambda_a-1)
                // m=2 for green strain case in Eq.(16)
                theta_ab=0.5*(epos[a]-epos[b])/(eigval[a]-eigval[b]);// Eq.(21)-1
            }
            ProjPos+=theta_ab*(Gab+Gba);
        }
    }
    return ProjPos;
}