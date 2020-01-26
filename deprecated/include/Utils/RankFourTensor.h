#ifndef ASFEM_RANKFOURTENSOR_H
#define ASFEM_RANKFOURTENSOR_H

#include <iostream>
#include <iomanip>
#include <cmath>


#include "MessagePrint/MessagePrint.h"

#include "RankTwoTensor.h"


using namespace std;

class RankTwoTensor;

class RankFourTensor
{
public:
    enum InitMethod{
        InitZero,
        InitIdentity,
        InitIdentity4,
        InitIdentitySymmetric4
    };
public:
    RankFourTensor(double val=0.0);
    RankFourTensor(InitMethod method=InitMethod::InitZero);

    inline int GetDim()const{return _N;}


    inline double GetIKjlComponent(const int &i,const int &k,
                                   const Eigen::Vector3d &grad_test,
                                   const Eigen::Vector3d &grad_phi)const{
        // function to get K_ik=C_ijkl*N,j*N,l(for assemble local K use rank-4 tensor!!!)
        // where j is the index of trial fun
        //       i is the index of test fun
        return ((*this)(i,1,k,1)*grad_phi.coeff(0)
               +(*this)(i,1,k,2)*grad_phi.coeff(1)
               +(*this)(i,1,k,3)*grad_phi.coeff(2))*grad_test.coeff(0)
              +((*this)(i,2,k,1)*grad_phi.coeff(0)
               +(*this)(i,2,k,2)*grad_phi.coeff(1)
               +(*this)(i,2,k,3)*grad_phi.coeff(2))*grad_test.coeff(1)
              +((*this)(i,3,k,1)*grad_phi.coeff(0)
               +(*this)(i,3,k,2)*grad_phi.coeff(1)
               +(*this)(i,3,k,3)*grad_phi.coeff(2))*grad_test.coeff(2);
    }

    //********************************************
    //*** For operator overload
    //********************************************
    inline double operator()(const int &i,const int &j,const int &k,const int &l) const{
        return _vals[(((i-1)*_N+j-1)*_N+k-1)*_N+l-1];
    }
    inline double& operator()(const int &i,const int &j,const int &k,const int &l){
        return _vals[(((i-1)*_N+j-1)*_N+k-1)*_N+l-1];
    }
    inline double  operator[](const int &i) const{
        return _vals[i-1];
    }
    inline double& operator[](const int &i){
        return _vals[i-1];
    }
    //************
    //*** =
    //************
    inline RankFourTensor& operator=(const double &a){
        for(int i=0;i<_N4;++i) _vals[i]=a;
        return *this;
    }
    inline RankFourTensor& operator=(const RankFourTensor &a){
        for(int i=0;i<_N4;++i) _vals[i]=a._vals[i];
        return *this;
    }
    //************
    //*** +
    //************
    inline RankFourTensor operator+(const double &a) const{
        RankFourTensor temp(0.0);
        for(int i=0;i<_N4;++i) temp._vals[i]=_vals[i]+a;
        return temp;
    }
    inline RankFourTensor operator+(const RankFourTensor &a) const{
        RankFourTensor temp(0.0);
        for(int i=0;i<_N4;++i) temp._vals[i]=_vals[i]+a._vals[i];
        return temp;
    }
    //**** for +=
    inline RankFourTensor& operator+=(const double &a){
        for(int i=0;i<_N4;++i) _vals[i]=_vals[i]+a;
        return *this;
    }
    inline RankFourTensor& operator+=(const RankFourTensor &a){
        for(int i=0;i<_N4;++i) _vals[i]=_vals[i]+a._vals[i];
        return *this;
    }
    //************
    //*** -
    //************
    inline RankFourTensor operator-(const double &a) const{
        RankFourTensor temp(0.0);
        for(int i=0;i<_N4;++i) temp._vals[i]=_vals[i]-a;
        return temp;
    }
    inline RankFourTensor operator-(const RankFourTensor &a) const{
        RankFourTensor temp(0.0);
        for(int i=0;i<_N4;++i) temp._vals[i]=_vals[i]-a._vals[i];
        return temp;
    }
    //**** for +=
    inline RankFourTensor& operator-=(const double &a){
        for(int i=0;i<_N4;++i) _vals[i]=_vals[i]-a;
        return *this;
    }
    inline RankFourTensor& operator-=(const RankFourTensor &a){
        for(int i=0;i<_N4;++i) _vals[i]=_vals[i]-a._vals[i];
        return *this;
    }
    //************
    //*** *
    //************
    inline RankFourTensor operator*(const double &a) const{
        RankFourTensor temp(0.0);
        for(int i=0;i<_N4;++i) temp._vals[i]=_vals[i]*a;
        return temp;
    }
    //*** be careful:
    // this must be put into cpp file, otherwise the compiler will complain the errors!!!
    friend RankFourTensor operator*(const double &lhs,const RankFourTensor &a);
    friend RankTwoTensor operator*(const RankTwoTensor &lhs,const RankFourTensor &a);
    //*** for mixed case, please put it into cpp
    RankTwoTensor operator*(const RankTwoTensor &a) const;

    inline RankFourTensor operator*(const RankFourTensor &a) const{
        // C_ijkl=A_ijmn*B_mnkl
        RankFourTensor temp(0.0);
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                for(int k=1;k<=_N;++k){
                    for(int l=1;l<=_N;++l){
                        temp(i,j,k,l)=0.0;
                        for(int m=1;m<=_N;++m){
                            for(int n=1;n<=_N;++n){
                                temp(i,j,k,l)+=(*this)(i,j,m,n)*a(m,n,k,l);
                            }
                        }
                    }
                }
            }
        }
        return temp;
    }
    //**** for *=
    inline RankFourTensor& operator*=(const double &a){
        for(int i=0;i<_N4;++i) _vals[i]=_vals[i]*a;
        return *this;
    }
    inline RankFourTensor& operator*=(const RankFourTensor &a){
        RankFourTensor temp(0.0);
        temp=(*this)*a;
        (*this)=temp;
        return *this;
    }
    //**********************************************
    //*** some setting functions
    //**********************************************
    inline void SetToZeros(){
        for(int i=0;i<_N4;++i) _vals[i]=0.0;
    }
    inline void SetToIdentity(){
        SetToZeros();
        for(int i=1;i<=_N;++i) (*this)(i,i,i,i)=1.0;
    }
    inline void SetToIdentity4(){
        SetToZeros();
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                for(int k=1;k<=_N;++k){
                    for(int l=1;l<=_N;++l){
                        (*this)(i,j,k,l)=1.0*((i==k)*(j==l));
                    }
                }
            }
        }
    }
    inline void SetToIdentitySymmetric4(){
        SetToZeros();
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                for(int k=1;k<=_N;++k){
                    for(int l=1;l<=_N;++l){
                        (*this)(i,j,k,l)=0.5*((i==k)&&(j==l))
                                        +0.5*((i==l)&&(j==k));
                    }
                }
            }
        }
    }
    //**********************************************
    //*** some fill-in functions
    //**********************************************
    inline int Eps(const int &i,const int &j) const{
        if(i==1&&j==2) {
            return 1;
        }
        else if(i==2&&j==1){
            return -1;
        }
        return 0;
    }
    inline int Eps(const int &i,const int &j,const int &k) const{
        if(i==1&&j>1&&k>1){
            return Eps(j-1,k-1);
        }
        else if(j==1&&i>1&&k>1){
            return -Eps(i-1,k-1);
        }
        else if(k==1&&i>1&&j>1){
            return Eps(i-1,j-1);
        }
        return 0;
    }
    //*** fill method
    void SetFromLameandG(const double &Lame,const double &G);
    void SetFromEandNu(const double &E,const double &Nu);
    void SetFromKandG(const double &K,const double &G);
    void SetFromSymmetric9(const vector<double> &vec);


    //*****************************************
    //*** Print rank-4 tensor
    //*****************************************
    inline void Print() const{
        for(int i=1;i<=_N;++i){
            for(int j=1;j<=_N;++j){
                printf("*** i=%d,j=%d: %14.6e  %14.6e  %14.6e***\n",i,j,(*this)(i,j,1,1),(*this)(i,j,1,2),(*this)(i,j,1,3));
                printf("***          %14.6e  %14.6e  %14.6e***\n",(*this)(i,j,2,1),(*this)(i,j,2,2),(*this)(i,j,2,3));
                printf("***          %14.6e  %14.6e  %14.6e***\n",(*this)(i,j,3,1),(*this)(i,j,3,2),(*this)(i,j,3,3));
            }
        }
    }
    //*******************
    inline void PrintVoigt() const{
        // taken from: https://en.wikipedia.org/wiki/Linear_elasticity
        printf("*** %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e***\n",(*this)(1,1,1,1),(*this)(1,1,2,2),(*this)(1,1,3,3),(*this)(1,1,2,3),(*this)(1,1,1,3),(*this)(1,1,1,2));
        printf("*** %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e***\n",(*this)(2,2,1,1),(*this)(2,2,2,2),(*this)(2,2,3,3),(*this)(2,2,2,3),(*this)(2,2,1,3),(*this)(2,2,1,2));
        printf("*** %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e***\n",(*this)(3,3,1,1),(*this)(3,3,2,2),(*this)(3,3,3,3),(*this)(3,3,2,3),(*this)(3,3,1,3),(*this)(3,3,1,2));
        printf("*** %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e***\n",(*this)(2,3,1,1),(*this)(2,3,2,2),(*this)(2,3,3,3),(*this)(2,3,2,3),(*this)(2,3,1,3),(*this)(2,3,1,2));
        printf("*** %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e***\n",(*this)(3,1,1,1),(*this)(3,1,2,2),(*this)(3,1,3,3),(*this)(3,1,2,3),(*this)(3,1,1,3),(*this)(3,1,1,2));
        printf("*** %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e***\n",(*this)(1,2,1,1),(*this)(1,2,2,2),(*this)(1,2,3,3),(*this)(1,2,2,3),(*this)(1,2,1,3),(*this)(1,2,1,2));
    }
private:
    double _vals[81];
    const int _N=3,_N2=9,_N4=81;
};

#endif //ASFEM_RANKFOURTENSOR_H