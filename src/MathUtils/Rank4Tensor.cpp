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
//+++ Purpose: Implement rank-4 tensor class for the common
//+++          tensor manipulation in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MathUtils/Rank4Tensor.h"

Rank4Tensor::Rank4Tensor(){
    m_vals.resize(81,0.0);
}
Rank4Tensor::Rank4Tensor(const double &val){
    m_vals.resize(81,val);
}
Rank4Tensor::Rank4Tensor(const Rank4Tensor &a){
    m_vals.resize(81,0.0);
    for(int i=0;i<N4;i++) m_vals[i]=a.m_vals[i];
}
Rank4Tensor::Rank4Tensor(const InitMethod &method){
    if(method==InitMethod::ZERO){
        m_vals.resize(81,0.0);
    }
    else if(method==InitMethod::IDENTITY){
        m_vals.resize(81,0.0);
        setToIdentity();
    }
    else if(method==InitMethod::IDENTITY4){
        m_vals.resize(81,0.0);
        setToIdentity4();
    }
    else if(method==InitMethod::IDENTITY4SYMMETRIC){
        m_vals.resize(81,0.0);
        setToIdentity4Symmetric();
    }
    else if(method==InitMethod::IDENTITY4TRANS){
        m_vals.resize(81,0.0);
        setIdentity4Transpose();
    }
    else{
        MessagePrinter::printErrorTxt("unsupported initialize method in rank-4 tensor");
        MessagePrinter::exitAsFem();
    }
}
Rank4Tensor::~Rank4Tensor(){
    m_vals.clear();
}
//************************************************************************
double Rank4Tensor::getVoigtComponent(const int &i,const int &j)const{
    if(i==1){
        if(j==1){
            return (*this)(1,1,1,1);
        }
        else if(j==2){
            return (*this)(1,1,2,2);
        }
        else if(j==3){
            return (*this)(1,1,3,3);
        }
        else if(j==4){
            return (*this)(1,1,2,3);
        }
        else if(j==5){
            return (*this)(1,1,1,3);
        }
        else if(j==6){
            return (*this)(1,1,1,2);
        }
        else{
            MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
    }
    else if(i==2){
        if(j==1){
            return (*this)(2,2,1,1);
        }
        else if(j==2){
            return (*this)(2,2,2,2);
        }
        else if(j==3){
            return (*this)(2,2,3,3);
        }
        else if(j==4){
            return (*this)(2,2,2,3);
        }
        else if(j==5){
            return (*this)(2,2,1,3);
        }
        else if(j==6){
            return (*this)(2,2,1,2);
        }
        else{
            MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
    }
    else if(i==3){
        if(j==1){
            return (*this)(3,3,1,1);
        }
        else if(j==2){
            return (*this)(3,3,2,2);
        }
        else if(j==3){
            return (*this)(3,3,3,3);
        }
        else if(j==4){
            return (*this)(3,3,2,3);
        }
        else if(j==5){
            return (*this)(3,3,1,3);
        }
        else if(j==6){
            return (*this)(3,3,1,2);
        }
        else{
            MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
    }
    else if(i==4){
        if(j==1){
            return (*this)(2,3,1,1);
        }
        else if(j==2){
            return (*this)(2,3,2,2);
        }
        else if(j==3){
            return (*this)(2,3,3,3);
        }
        else if(j==4){
            return (*this)(2,3,2,3);
        }
        else if(j==5){
            return (*this)(2,3,1,3);
        }
        else if(j==6){
            return (*this)(2,3,1,2);
        }
        else{
            MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
    }
    else if(i==5){
        if(j==1){
            return (*this)(1,3,1,1);
        }
        else if(j==2){
            return (*this)(1,3,2,2);
        }
        else if(j==3){
            return (*this)(1,3,3,3);
        }
        else if(j==4){
            return (*this)(1,3,2,3);
        }
        else if(j==5){
            return (*this)(1,3,1,3);
        }
        else if(j==6){
            return (*this)(1,3,1,2);
        }
        else{
            MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
    }
    else if(i==6){
        if(j==1){
            return (*this)(1,2,1,1);
        }
        else if(j==2){
            return (*this)(1,2,2,2);
        }
        else if(j==3){
            return (*this)(1,2,3,3);
        }
        else if(j==4){
            return (*this)(1,2,2,3);
        }
        else if(j==5){
            return (*this)(1,2,1,3);
        }
        else if(j==6){
            return (*this)(1,2,1,2);
        }
        else{
            MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
        MessagePrinter::exitAsFem();
    }
    return 0;
}
//************************************************************************
double& Rank4Tensor::voigtComponent(const int &i,const int &j){
    if(i==1){
        if(j==1){
            return (*this)(1,1,1,1);
        }
        else if(j==2){
            return (*this)(1,1,2,2);
        }
        else if(j==3){
            return (*this)(1,1,3,3);
        }
        else if(j==4){
            return (*this)(1,1,2,3);
        }
        else if(j==5){
            return (*this)(1,1,1,3);
        }
        else if(j==6){
            return (*this)(1,1,1,2);
        }
        else{
            MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
    }
    else if(i==2){
        if(j==1){
            return (*this)(2,2,1,1);
        }
        else if(j==2){
            return (*this)(2,2,2,2);
        }
        else if(j==3){
            return (*this)(2,2,3,3);
        }
        else if(j==4){
            return (*this)(2,2,2,3);
        }
        else if(j==5){
            return (*this)(2,2,1,3);
        }
        else if(j==6){
            return (*this)(2,2,1,2);
        }
        else{
            MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
    }
    else if(i==3){
        if(j==1){
            return (*this)(3,3,1,1);
        }
        else if(j==2){
            return (*this)(3,3,2,2);
        }
        else if(j==3){
            return (*this)(3,3,3,3);
        }
        else if(j==4){
            return (*this)(3,3,2,3);
        }
        else if(j==5){
            return (*this)(3,3,1,3);
        }
        else if(j==6){
            return (*this)(3,3,1,2);
        }
        else{
            MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
    }
    else if(i==4){
        if(j==1){
            return (*this)(2,3,1,1);
        }
        else if(j==2){
            return (*this)(2,3,2,2);
        }
        else if(j==3){
            return (*this)(2,3,3,3);
        }
        else if(j==4){
            return (*this)(2,3,2,3);
        }
        else if(j==5){
            return (*this)(2,3,1,3);
        }
        else if(j==6){
            return (*this)(2,3,1,2);
        }
        else{
            MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
    }
    else if(i==5){
        if(j==1){
            return (*this)(1,3,1,1);
        }
        else if(j==2){
            return (*this)(1,3,2,2);
        }
        else if(j==3){
            return (*this)(1,3,3,3);
        }
        else if(j==4){
            return (*this)(1,3,2,3);
        }
        else if(j==5){
            return (*this)(1,3,1,3);
        }
        else if(j==6){
            return (*this)(1,3,1,2);
        }
        else{
            MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
    }
    else if(i==6){
        if(j==1){
            return (*this)(1,2,1,1);
        }
        else if(j==2){
            return (*this)(1,2,2,2);
        }
        else if(j==3){
            return (*this)(1,2,3,3);
        }
        else if(j==4){
            return (*this)(1,2,2,3);
        }
        else if(j==5){
            return (*this)(1,2,1,3);
        }
        else if(j==6){
            return (*this)(1,2,1,2);
        }
        else{
            MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        MessagePrinter::printErrorTxt("i="+to_string(i)+", j="+to_string(j)+" is invalid for the voigt notation in rank-4 tensor");
        MessagePrinter::exitAsFem();
    }
    return (*this)(1,1,1,1);
}
//************************************************************************
Rank4Tensor Rank4Tensor::operator*(const Rank2Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=0.0;
                    for(int m=1;m<=N;m++){
                        temp(i,j,k,l)+=(*this)(i,j,k,m)*a(m,l);
                    }
                }
            }
        }
    }
    return temp;
}
Rank2Tensor Rank4Tensor::doubledot(const Rank2Tensor &a) const{
    Rank2Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            temp(i,j)=0.0;
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j)+=(*this)(i,j,k,l)*a(k,l);
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank4Tensor::doubledot(const Rank4Tensor &a) const{
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=0.0;
                    for(int m=1;m<=N;m++){
                        for(int n=1;n<=N;n++){
                            temp(i,j,k,l)+=(*this)(i,j,m,n)*a(m,n,k,l);
                        }
                    }
                }
            }
        }
    }
    return temp;
}
//*******************************************************************
//*** for left hand side manipulation
//*******************************************************************
Rank4Tensor operator*(const double &lhs,const Rank4Tensor &a){
    Rank4Tensor temp(0.0);
    for(int i=0;i<a.N4;i++) temp.m_vals[i]=lhs*a.m_vals[i];
    return temp;
}
Rank4Tensor operator*(const Rank2Tensor &lhs,const Rank4Tensor &a){
    // C_ijkl=B_ip*A_pjkl
    Rank4Tensor temp(0.0);
    for(int i=1;i<=a.N;i++){
        for(int j=1;j<=a.N;j++){
            for(int k=1;k<=a.N;k++){
                for(int l=1;l<=a.N;l++){
                    temp(i,j,k,l)=0.0;
                    for(int p=1;p<=a.N;p++){
                        temp(i,j,k,l)+=lhs(i,p)*a(p,j,k,l);
                    }
                }
            }
        }
    }
    return temp;
}
//*******************************************************************
//*** fill up the rank-4 tensor
//*******************************************************************
void Rank4Tensor::setFromLameAndG(const double &Lame,const double &G){
    // taken from: https://en.wikipedia.org/wiki/Linear_elasticity
    // C_ijkl = Lame*de_ij*de_kl + G*(de_ik*de_jl + de_il*de_jk)
    setToZeros();
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    (*this)(i,j,k,l)=Lame*(i==j)*(k==l)
                                    +G*(i==k)*(j==l)
                                    +G*(i==l)*(j==k);
                }
            }
        }
    }
}
void Rank4Tensor::setFromEAndNu(const double &E,const double &Nu){
    // taken from:  https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
    double Lame=E*Nu/((1.0+Nu)*(1.0-2.0*Nu));
    double G=E/(2.0*(1.0+Nu));
    setFromLameAndG(Lame,G);
}
void Rank4Tensor::setFromKAndG(const double &K,const double &G){
    // taken from:  https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
    double Lame=K-2.0*G/3.0;
    setFromLameAndG(Lame,G);
}
void Rank4Tensor::setFromSymmetric9(const vector<double> &vec){
    if(vec.size()<9){
        MessagePrinter::printErrorTxt("Symmetric9 fill method for rank-4 tensor need at least 9 elements!!!");
        MessagePrinter::exitAsFem();
    }
    //C1111  C1122  C1133   0     0     0
    // 0     C2222  C2233   0     0     0
    // 0      0     C3333   0     0     0
    // 0      0      0     C2323  0     0
    // 0      0      0      0    C1313  0
    // 0      0      0      0     0    C1212

    // C1111,C1122,C1133,C2222,C2233,C3333,C2323,C1313,C1212
    // C11,  C12,  C13,  C22,  C23,  C33,  C44,  C55,  C66
    setToZeros();
    (*this)(1,1,1,1)=vec[0];
    (*this)(1,1,2,2)=vec[1];
    (*this)(1,1,3,3)=vec[2];
    (*this)(2,2,2,2)=vec[3];
    (*this)(2,2,3,3)=vec[4];
    (*this)(3,3,3,3)=vec[5];
    (*this)(2,3,2,3)=vec[6];
    (*this)(1,3,1,3)=vec[7];
    (*this)(1,2,1,2)=vec[8];

    // complete the symmetric part
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    (*this)(i,j,l,k)=(*this)(i,j,k,l);

                    (*this)(j,i,k,l)=(*this)(i,j,k,l);
                    (*this)(j,i,l,k)=(*this)(i,j,k,l);

                    (*this)(k,l,i,j)=(*this)(i,j,k,l);
                    (*this)(k,l,j,i)=(*this)(i,j,k,l);

                    (*this)(l,k,j,i)=(*this)(i,j,k,l);
                    (*this)(l,k,i,j)=(*this)(i,j,k,l);
                }
            }
        }
    }
}
void Rank4Tensor::setToOrthotropic(const vector<double> &vec){
    if(vec.size()<9){
        MessagePrinter::printErrorTxt("Orthotropic fill method for rank-4 tensor need at least 9 elements!!!");
        MessagePrinter::exitAsFem();
    }
    //C1111  C1122  C1133   0     0     0
    // 0     C2222  C2233   0     0     0
    // 0      0     C3333   0     0     0
    // 0      0      0     C2323  0     0
    // 0      0      0      0    C1313  0
    // 0      0      0      0     0    C1212

    // C1111,C1122,C1133,C2222,C2233,C3333,C2323,C1313,C1212
    // C11,  C12,  C13,  C22,  C23,  C33,  C44,  C55,  C66
    setToZeros();
    (*this)(1,1,1,1)=vec[0];
    (*this)(1,1,2,2)=vec[1];
    (*this)(1,1,3,3)=vec[2];
    (*this)(2,2,2,2)=vec[3];
    (*this)(2,2,3,3)=vec[4];
    (*this)(3,3,3,3)=vec[5];
    (*this)(2,3,2,3)=vec[6];
    (*this)(1,3,1,3)=vec[7];
    (*this)(1,2,1,2)=vec[8];

    // complete the symmetric part
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    (*this)(i,j,l,k)=(*this)(i,j,k,l);

                    (*this)(j,i,k,l)=(*this)(i,j,k,l);
                    (*this)(j,i,l,k)=(*this)(i,j,k,l);

                    (*this)(k,l,i,j)=(*this)(i,j,k,l);
                    (*this)(k,l,j,i)=(*this)(i,j,k,l);

                    (*this)(l,k,j,i)=(*this)(i,j,k,l);
                    (*this)(l,k,i,j)=(*this)(i,j,k,l);
                }
            }
        }
    }
}
//**********************************************************************
//*** for advanced mathematic manipulation
//**********************************************************************
Rank4Tensor Rank4Tensor::rotate(const Rank2Tensor &rot) const{
    //C_ijkl=C_mnpq*R_im*R_jn*R_kp*R_lq
    Rank4Tensor temp(0.0);
    for(int i=1;i<=3;i++){
        for(int j=1;j<=3;j++){
            for(int k=1;k<=3;k++){
                for(int l=1;l<=3;l++){
                    temp(i,j,k,l)=0.0;
                    for(int m=1;m<=3;m++){
                        for(int n=1;n<=3;n++){
                            for(int p=1;p<=3;p++){
                                for(int q=1;q<=3;q++){
                                    temp(i,j,k,l)+=(*this)(m,n,p,q)*rot(i,m)*rot(j,n)*rot(k,p)*rot(l,q);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return temp;
}
void Rank4Tensor::rotated(const Rank2Tensor &rot){
    Rank4Tensor temp(0.0);
    temp=(*this).rotate(rot);
    (*this)=temp;
}
Rank4Tensor Rank4Tensor::pushForward(const Rank2Tensor &F) const{
    // push forward the jacobian from reference one to the current configuration or similar operation
    // new Ran-4 tensor c_ijkl=F_iA*F_jB*F_kC*F_lD*C_ABCD
    Rank4Tensor temp(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    temp(i,j,k,l)=0.0;
                    for(int A=1;A<=N;A++){
                        for(int B=1;B<=N;B++){
                            for(int C=1;C<=N;C++){
                                for(int D=1;D<=N;D++){
                                    temp(i,j,k,l)+=F(i,A)*F(j,B)*F(k,C)*F(l,D)*(*this)(A,B,C,D);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return temp;
}
Rank4Tensor Rank4Tensor::conjPushForward(const Rank2Tensor &F) const{
    // calculate the conjugate rank-4 tensor regulate by F
    // C_ijkl=F_im*C_mjnl*F_kn
    Rank4Tensor A(0.0);
    for(int i=1;i<=N;i++){
        for(int j=1;j<=N;j++){
            for(int k=1;k<=N;k++){
                for(int l=1;l<=N;l++){
                    for(int m=1;m<=N;m++){
                        for(int n=1;n<=N;n++){
                            A(i,j,k,l)+=F(i,m)*(*this)(m,j,n,l)*F(k,n);
                        }
                    }
                }
            }
        }
    }

    return A;
}