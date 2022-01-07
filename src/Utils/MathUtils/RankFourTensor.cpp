//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.10.17
//+++ Purpose: Implement rank-4 tensor class for some common
//+++          tensor calculation for solid-mechanics problem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Utils/RankFourTensor.h"


// for the constructor in different purpose
RankFourTensor::RankFourTensor()
:_N(3),_N2(3*3),_N4(3*3*3*3){
    for(int i=0;i<_N4;i++) _vals[i]=0.0;
}
RankFourTensor::RankFourTensor(const double &val)
:_N(3),_N2(3*3),_N4(3*3*3*3){
    for(int i=0;i<_N4;i++) _vals[i]=val;
}
RankFourTensor::RankFourTensor(const RankFourTensor &a)
:_N(3),_N2(3*3),_N4(3*3*3*3){
    for(int i=0;i<_N4;i++) _vals[i]=a._vals[i];
}
//*********************
RankFourTensor::RankFourTensor(const InitMethod &method)
:_N(3),_N2(3*3),_N4(3*3*3*3){
    switch (method)
    {
    case InitMethod::InitZero:
        SetToZeros();
        break;
    case InitMethod::InitIdentity:
        SetToIdentity();
        break;
    case InitMethod::InitIdentity4:
        SetToIdentity4();
        break;
    case InitMethod::InitIdentitySymmetric4:
        SetToIdentitySymmetric4();
        break;
    case InitMethod::InitRandom:
        SetToRandom();
        break;
    default:
        MessagePrinter::PrintErrorTxt("unsupported initializing method for rank-4 tensor");
        MessagePrinter::AsFem_Exit();
        break;
    }
}

double RankFourTensor::VoigtIJcomponent(const int &i,const int &j) const{
    if(i<1||i>6 || j<1||j>6){
        MessagePrinter::PrintErrorTxt("your i or j (1~6) is out of range when you call a rank-4 tensor in voigt notation");
        MessagePrinter::AsFem_Exit();
    }
    if(i==1 && j==1){
        return (*this)(1,1,1,1);
    }
    else if(i==1 && j==2){
        return (*this)(1,1,2,2);
    }
    else if(i==1 && j==3){
        return (*this)(1,1,3,3);
    }
    else if(i==1 && j==4){
        return (*this)(1,1,2,3);
    }
    else if(i==1 && j==5){
        return (*this)(1,1,1,3);
    }
    else if(i==1 && j==6){
        return (*this)(1,1,1,2);
    }
    //  
    else if(i==2 && j==1){
        return (*this)(2,2,1,1);
    }
    else if(i==2 && j==2){
        return (*this)(2,2,2,2);
    }
    else if(i==2 && j==3){
        return (*this)(2,2,3,3);
    }
    else if(i==2 && j==4){
        return (*this)(2,2,2,3);
    }
    else if(i==2 && j==5){
        return (*this)(2,2,1,3);
    }
    else if(i==2 && j==6){
        return (*this)(2,2,1,2);
    }
    //  
    else if(i==3 && j==1){
        return (*this)(3,3,1,1);
    }
    else if(i==3 && j==2){
        return (*this)(3,3,2,2);
    }
    else if(i==3 && j==3){
        return (*this)(3,3,3,3);
    }
    else if(i==3 && j==4){
        return (*this)(3,3,2,3);
    }
    else if(i==3 && j==5){
        return (*this)(3,3,1,3);
    }
    else if(i==3 && j==6){
        return (*this)(3,3,1,2);
    }
    //  
    else if(i==4 && j==1){
        return (*this)(2,3,1,1);
    }
    else if(i==4 && j==2){
        return (*this)(2,3,2,2);
    }
    else if(i==4 && j==3){
        return (*this)(2,3,3,3);
    }
    else if(i==4 && j==4){
        return (*this)(2,3,2,3);
    }
    else if(i==4 && j==5){
        return (*this)(2,3,1,3);
    }
    else if(i==4 && j==6){
        return (*this)(2,3,1,2);
    }
    //  
    else if(i==5 && j==1){
        return (*this)(3,1,1,1);
    }
    else if(i==5 && j==2){
        return (*this)(3,1,2,2);
    }
    else if(i==5 && j==3){
        return (*this)(3,1,3,3);
    }
    else if(i==5 && j==4){
        return (*this)(3,1,2,3);
    }
    else if(i==5 && j==5){
        return (*this)(3,1,1,3);
    }
    else if(i==5 && j==6){
        return (*this)(3,1,1,2);
    }
    //  
    else if(i==6 && j==1){
        return (*this)(1,2,1,1);
    }
    else if(i==6 && j==2){
        return (*this)(1,2,2,2);
    }
    else if(i==6 && j==3){
        return (*this)(1,2,3,3);
    }
    else if(i==6 && j==4){
        return (*this)(1,2,2,3);
    }
    else if(i==6 && j==5){
        return (*this)(1,2,1,3);
    }
    else if(i==6 && j==6){
        return (*this)(1,2,1,2);
    }
       
    return 0.0;
}
//****************************************
//*** for fill-in method
//****************************************
void RankFourTensor::SetFromLameandG(const double &Lame,const double &G){
    // taken from: https://en.wikipedia.org/wiki/Linear_elasticity
    // C_ijkl = Lame*de_ij*de_kl + G*(de_ik*de_jl + de_il*de_jk)
    SetToZeros();
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
                    (*this)(i,j,k,l)=Lame*(i==j)*(k==l)
                                    +G*(i==k)*(j==l)
                                    +G*(i==l)*(j==k);
                }
            }
        }
    }
}
//*** use Youngs modulus and poisson ratio
void RankFourTensor::SetFromEandNu(const double &E,const double &Nu){
    // taken from:  https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
    double Lame=E*Nu/((1.0+Nu)*(1.0-2.0*Nu));
    double G=E/(2.0*(1.0+Nu));
    SetFromLameandG(Lame,G);
}
void RankFourTensor::SetFromKandG(const double &K,const double &G){
    // taken from:  https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
    double Lame=K-2.0*G/3.0;
    SetFromLameandG(Lame,G);
}
//*** use 9-independent element to fill the rank-4 tensor
void RankFourTensor::SetFromSymmetric9(const vector<double> &vec){
    if(vec.size()<9){
        MessagePrinter::PrintErrorTxt("Symmetric9 fill method for rank-4 tensor need at least 9 elements!!!");
        MessagePrinter::AsFem_Exit();
    }
    //C1111  C1122  C1133   0     0     0
    // 0     C2222  C2233   0     0     0
    // 0      0     C3333   0     0     0
    // 0      0      0     C2323  0     0
    // 0      0      0      0    C1313  0
    // 0      0      0      0     0    C1212
    // C1111,C1122,C1133,C2222,C2233,C3333,C2323,C1313,C1212
    // C11,  C12,  C13,  C22,  C23,  C33,  C44,  C55,  C66
    SetToZeros();
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
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
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
void RankFourTensor::SetToOrthotropic(const vector<double> &vec){
    if(vec.size()<9){
        MessagePrinter::PrintErrorTxt("Orthotropic fill method for rank-4 tensor need at least 9 elements!!!");
        MessagePrinter::AsFem_Exit();
    }
    //C1111  C1122  C1133   0     0     0
    // 0     C2222  C2233   0     0     0
    // 0      0     C3333   0     0     0
    // 0      0      0     C2323  0     0
    // 0      0      0      0    C1313  0
    // 0      0      0      0     0    C1212
    // C1111,C1122,C1133,C2222,C2233,C3333,C2323,C1313,C1212
    // C11,  C12,  C13,  C22,  C23,  C33,  C44,  C55,  C66
    SetToZeros();
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
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
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
//*** for mixed case *
RankFourTensor RankFourTensor::operator*(const RankTwoTensor &a) const{
    RankFourTensor temp(0.0);
    for(int i=1;i<=_N;i++){
        for(int j=1;j<=_N;j++){
            for(int k=1;k<=_N;k++){
                for(int l=1;l<=_N;l++){
                    for(int m=1;m<=_N;m++){
                        temp(i,j,k,l)+=(*this)(i,j,k,m)*a(m,l);
                    }
                }
            }
        }
    }
    return temp;
}
//**** for left hand scale value time rank-4 tensor
RankFourTensor operator*(const double &lhs,const RankFourTensor &a) {
    RankFourTensor temp(0.0);
    for(int i=0;i<a._N4;i++) temp._vals[i]=lhs*a._vals[i];
    return temp;
}
//**** for left hand rank-2 tensor times rank-4 tensor(still return rank-2 tensor)
RankFourTensor operator*(const RankTwoTensor &lhs,const RankFourTensor &a){
    // C_ijkl=B_ip*A_pjkl
    RankFourTensor temp(0.0);
    for(int i=1;i<=a._N;++i){
        for(int j=1;j<=a._N;++j){
            for(int k=1;k<=a._N;++k){
                for(int l=1;l<=a._N;++l){
                    for(int p=1;p<=a._N;++p){
                        temp(i,j,k,l)+=lhs(i,p)*a(p,j,k,l);
                    }
                }
            }
        }
    }
    return temp;
}
// for double dot operator
RankTwoTensor RankFourTensor::DoubleDot(const RankTwoTensor &a) const{
    // A_ijkl:B_kl = Cij
    RankTwoTensor temp(0.0);
    for(int i=1;i<=_N;i++){
        for(int j=1;j<=_N;j++){
            for(int k=1;k<=_N;k++){
                for(int l=1;l<=_N;l++){
                    temp(i,j)+=(*this)(i,j,k,l)*a(k,l);
                }
            }
        }
    }
    return temp;
}
//**************************************************
//*** For rotation of a rank-4 tensor by the rank-2
//*** rotation tensor
//**************************************************
RankFourTensor RankFourTensor::Rotate(const RankTwoTensor &rotate) const{
    //C_ijkl=C_mnpq*R_im*R_jn*R_kp*R_lq
    RankFourTensor temp(0.0);
    for(int i=1;i<=3;i++){
        for(int j=1;j<=3;j++){
            for(int k=1;k<=3;k++){
                for(int l=1;l<=3;l++){
                    for(int m=1;m<=3;m++){
                        for(int n=1;n<=3;n++){
                            for(int p=1;p<=3;p++){
                                for(int q=1;q<=3;q++){
                                    temp(i,j,k,l)+=(*this)(m,n,p,q)*rotate(i,m)*rotate(j,n)*rotate(k,p)*rotate(l,q);
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
//****************************************************
//*** For push forward operator
//****************************************************
RankFourTensor RankFourTensor::PushByF(const RankTwoTensor &F) const{
    // new Rank-4 tensor=F*R4Old*F^T
    RankFourTensor A(0.0),B(0.0);
    RankTwoTensor Ft(0.0);
    Ft=F.Transpose();

    A.SetToZeros();
    for(int i=1;i<=_N;i++){
        for(int j=1;j<=_N;j++){
            for(int k=1;k<=_N;k++){
                for(int l=1;l<=_N;l++){
                    for(int m=1;m<=_N;m++){
                        A(i,j,k,l)+=F(i,m)*(*this)(m,j,k,l);
                    }
                }
            }
        }
    }

    B.SetToZeros();
    for(int i=1;i<=_N;i++){
        for(int j=1;j<=_N;j++){
            for(int k=1;k<=_N;k++){
                for(int l=1;l<=_N;l++){
                    for(int m=1;m<=_N;m++){
                        B(i,j,k,l)+=A(i,j,k,m)*Ft(m,l);
                    }
                }
            }
        }
    }
    return B;
}
//*** for the push-forward from reference to current
RankFourTensor RankFourTensor::PushForward(const RankTwoTensor &F) const{
    // push forward the jacobian from reference one to the current configuration or similar operation
    // new Ran-4 tensor ijkl=F_iA*F_jB*F_kC*F_lD*C_ABCD
    RankFourTensor temp;
    temp.SetToZeros();
    for(int i=1;i<=_N;i++){
        for(int j=1;j<=_N;j++){
            for(int k=1;k<=_N;k++){
                for(int l=1;l<=_N;l++){
                    for(int A=1;A<=_N;A++){
                        for(int B=1;B<=_N;B++){
                            for(int C=1;C<=_N;C++){
                                for(int D=1;D<=_N;D++){
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
