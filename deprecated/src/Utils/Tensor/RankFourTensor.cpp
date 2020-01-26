#include "Utils/RankFourTensor.h"

RankFourTensor::RankFourTensor(double val)
:_N(3),_N2(3*3),_N4(3*3*3*3){
    for(int i=0;i<_N4;++i) _vals[i]=val;
}
//*********************
RankFourTensor::RankFourTensor(InitMethod method)
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
    default:
        cout<<"*** Error: unsupported init method for rank-4 tensor !!!   ***"<<endl;
        Msg_AsFem_Exit();
        break;
    }
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
        cout<<"*** Error: Symmetric9 rank-9 tensor need 9 elements!!!     ***"<<endl;
        Msg_AsFem_Exit();
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
//*************************************
//*** for Rank-4*Rank-2 calculation ***
//*************************************
RankTwoTensor RankFourTensor::operator*(const RankTwoTensor &a) const{
    // C_ijkl*a_kl
    RankTwoTensor temp(0.0);
    for(int i=1;i<=_N;i++){
        for(int j=1;j<=_N;++j){
            temp(i,j)=0.0;
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
                    temp(i,j)+=(*this)(i,j,k,l)*a(k,l);
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
RankTwoTensor operator*(const RankTwoTensor &lhs,const RankFourTensor &a){
    // b_ij=a_kl*c_klij
    RankTwoTensor b(0.0);
    for(int i=1;i<=a._N;++i){
        for(int j=1;j<=a._N;++j){
            b(i,j)=0.0;
            for(int k=1;k<=a._N;++k){
                for(int l=1;l<=a._N;++l){
                    b(i,j)+=lhs(k,l)*a(k,l,i,j);
                }
            }
        }
    }
    return b;
}