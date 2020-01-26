#include "Utils/RankTwoTensor.h"

//***************************************
//*** the constructors for different usage!!!
//***************************************
RankTwoTensor::RankTwoTensor(double val)
:_N(3),_N2(3*3){
    for(int i=0;i<_N2;++i){
        _vals[i]=val;
    }
}
//***
RankTwoTensor::RankTwoTensor(InitMethod method)
:_N(3),_N2(3*3){
    switch(method){
        case InitMethod::InitZero:
            SetToZeros();
            break;
        case InitMethod::InitIdentity:
            SetToIdentity();
            break;
        default:
            cout<<"*** Error: unsupported rank-2 tensor fill method !!!       ***"<<endl;
            Msg_AsFem_Exit();
            break;
    }
}
//***(for deformation gradient case)
RankTwoTensor::RankTwoTensor(const Eigen::Vector3d &r1,
                             const Eigen::Vector3d &r2)
:_N(3),_N2(3*3){
    (*this)(1,1)=r1.coeff(0);(*this)(1,2)=r1.coeff(1);(*this)(1,3)=r1.coeff(2);
    (*this)(2,1)=r2.coeff(0);(*this)(2,2)=r2.coeff(1);(*this)(2,3)=r2.coeff(2);

    (*this)(3,1)=0.0;(*this)(3,2)=0.0;(*this)(3,3)=0.0;
}
RankTwoTensor::RankTwoTensor(const Eigen::Vector3d &r1,
                             const Eigen::Vector3d &r2,
                             const Eigen::Vector3d &r3)
:_N(3),_N2(3*3){
    (*this)(1,1)=r1.coeff(0);(*this)(1,2)=r1.coeff(1);(*this)(1,3)=r1.coeff(2);
    (*this)(2,1)=r2.coeff(0);(*this)(2,2)=r2.coeff(1);(*this)(2,3)=r2.coeff(2);
    (*this)(3,1)=r3.coeff(0);(*this)(3,2)=r3.coeff(1);(*this)(3,3)=r3.coeff(2);
}
//****(from voigt notation)
RankTwoTensor::RankTwoTensor(const double &v11,const double &v22,const double &v12)
:_N(3),_N2(3*3){
    (*this)(1,1)=v11;(*this)(1,2)=v12;(*this)(1,3)=0.0;
    (*this)(2,1)=v12;(*this)(2,2)=v22;(*this)(2,3)=0.0;
    (*this)(3,1)=0.0;(*this)(3,2)=0.0;(*this)(3,3)=0.0;
}
RankTwoTensor::RankTwoTensor(const double &v11,const double &v22,const double &v33,
                             const double &v23,const double &v31,const double &v12)
:_N(3),_N2(3*3){
    (*this)(1,1)=v11;(*this)(1,2)=v12;(*this)(1,3)=v31;
    (*this)(2,1)=v12;(*this)(2,2)=v22;(*this)(2,3)=v23;
    (*this)(3,1)=v31;(*this)(3,2)=v23;(*this)(3,3)=v33;
}
RankTwoTensor::RankTwoTensor(const double &v11,const double &v12,
                             const double &v21,const double &v22)
:_N(3),_N2(3*3){
    // for 2d voigt
    (*this)(1,1)=v11;(*this)(1,2)=v12;(*this)(1,3)=0.0;
    (*this)(2,1)=v21;(*this)(2,2)=v22;(*this)(2,3)=0.0;
    (*this)(3,1)=0.0;(*this)(3,2)=0.0;(*this)(3,3)=0.0;
}
RankTwoTensor::RankTwoTensor(const double &v11,const double &v12,const double &v13,
                  const double &v21,const double &v22,const double &v23,
                  const double &v31,const double &v32,const double &v33)
:_N(3),_N2(3*3){
    // for 3d voigt
    (*this)(1,1)=v11;(*this)(1,2)=v12;(*this)(1,3)=v13;
    (*this)(2,1)=v21;(*this)(2,2)=v22;(*this)(2,3)=v23;
    (*this)(3,1)=v31;(*this)(3,2)=v32;(*this)(3,3)=v33;
}
//**********************************************************
//*** Some setting functions
//**********************************************************
void RankTwoTensor::SetFromGradU(const Eigen::Vector3d &gradUx,const Eigen::Vector3d &gradUy){
    (*this)(1,1)=gradUx.coeff(0);(*this)(1,2)=gradUx.coeff(1);(*this)(1,3)=gradUx.coeff(2);
    (*this)(2,1)=gradUy.coeff(0);(*this)(2,2)=gradUy.coeff(1);(*this)(2,3)=gradUy.coeff(2);
    (*this)(3,1)=            0.0;(*this)(3,2)=            0.0;(*this)(3,3)=0.0;
}
void RankTwoTensor::SetFromGradU(const Eigen::Vector3d &gradUx,
                                 const Eigen::Vector3d &gradUy,
                                 const Eigen::Vector3d &gradUz){
    (*this)(1,1)=gradUx.coeff(0);(*this)(1,2)=gradUx.coeff(1);(*this)(1,3)=gradUx.coeff(2);
    (*this)(2,1)=gradUy.coeff(0);(*this)(2,2)=gradUy.coeff(1);(*this)(2,3)=gradUy.coeff(2);
    (*this)(3,1)=gradUz.coeff(0);(*this)(3,2)=gradUz.coeff(1);(*this)(3,3)=gradUz.coeff(2);
}
//********************************************************
//set tensor from voigt notation
//*****
void RankTwoTensor::SetFromVoigt(const double &v11,const double &v22,const double &v12){
    (*this)(1,1)=v11;(*this)(1,2)=v12;(*this)(1,3)=0.0;
    (*this)(2,1)=v12;(*this)(2,2)=v22;(*this)(2,3)=0.0;
    (*this)(3,1)=0.0;(*this)(3,2)=0.0;(*this)(3,3)=0.0;
}
void RankTwoTensor::SetFromVoigt(const double &v11,const double &v22,const double &v33,
                                 const double &v23,const double &v31,const double &v12){
    (*this)(1,1)=v11;(*this)(1,2)=v12;(*this)(1,3)=v31;
    (*this)(2,1)=v12;(*this)(2,2)=v22;(*this)(2,3)=v23;
    (*this)(3,1)=v31;(*this)(3,2)=v23;(*this)(3,3)=v33;
}
//****** for vector cross dot vector, which gives you a rank-2 tensor
// RankTwoTensor CrossDot(const Eigen::Vector3d &a,const Eigen::Vector3d &b){
//     RankTwoTensor temp(0.0);
//     for(int i=1;i<=_N;++i){
//         for(int j=1;j<=_N;++j){
//             temp(i,j)=a.coeff(i-1)*b.coeff(j-1);
//         }
//     }
//     return temp;
// }
//****** for left hand scale times rank-2 tensor
RankTwoTensor operator*(const double &lhs,const RankTwoTensor &a){
    RankTwoTensor temp(0.0);
    for(int i=0;i<a._N2;++i) temp._vals[i]=lhs*a._vals[i];
    return temp;
}
//**** for left hand vector times rank-2 tensor(return vector)
Eigen::Vector3d operator*(const Eigen::Vector3d &lhs,const RankTwoTensor &a){
    Eigen::Vector3d temp;
    for(int j=1;j<=a._N;++j){
        temp.coeffRef(j-1)=0.0;
        for(int i=1;i<=a._N;++i){
            temp.coeffRef(j)+=lhs.coeff(i-1)*a(i,j);
        }
    }
    return temp;
}
//*******************************************************************
//*** some higher order tensor calculation
//*******************************************************************
RankFourTensor RankTwoTensor::CrossDot(const RankTwoTensor &a) const{
    // return C_ijkl=a_ij*b_kl
    RankFourTensor temp(0.0);
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
                    temp(i,j,k,l)=(*this)(i,j)*a(k,l);
                }
            }
        }
    }
    return temp;
}
//*** for mixed case
RankFourTensor RankTwoTensor::IJlkDot(const RankTwoTensor &a) const{
    RankFourTensor temp(0.0);
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
                    temp(i,j,k,l)=(*this)(i,j)*a(l,k);
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For IkJl and IklJ
//********************************************************
RankFourTensor RankTwoTensor::IkJlDot(const RankTwoTensor &a) const{
    RankFourTensor temp(0.0);
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
                    temp(i,j,k,l)=(*this)(i,k)*a(j,l);
                }
            }
        }
    }
    return temp;
}
RankFourTensor RankTwoTensor::IklJDot(const RankTwoTensor &a) const{
    RankFourTensor temp(0.0);
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
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
RankFourTensor RankTwoTensor::IlJkDot(const RankTwoTensor &a) const{
    RankFourTensor temp(0.0);
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
                    temp(i,j,k,l)=(*this)(i,l)*a(j,k);
                }
            }
        }
    }
    return temp;
}
RankFourTensor RankTwoTensor::IlkJDot(const RankTwoTensor &a) const{
    RankFourTensor temp(0.0);
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
                    temp(i,j,k,l)=(*this)(i,l)*a(k,j);
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For JkIl and JklI
//********************************************************
RankFourTensor RankTwoTensor::JkIlDot(const RankTwoTensor &a) const{
    RankFourTensor temp(0.0);
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
                    temp(i,j,k,l)=(*this)(j,k)*a(i,l);
                }
            }
        }
    }
    return temp;
}
RankFourTensor RankTwoTensor::JklIDot(const RankTwoTensor &a) const{
    RankFourTensor temp(0.0);
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
                    temp(i,j,k,l)=(*this)(j,k)*a(l,i);
                }
            }
        }
    }
    return temp;
}
//********************************************************
//*** For JkIl and JklI
//********************************************************
RankFourTensor RankTwoTensor::klJIDot(const RankTwoTensor &a) const{
    RankFourTensor temp(0.0);
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
                    temp(i,j,k,l)=(*this)(k,l)*a(j,i);
                }
            }
        }
    }
    return temp;
}
RankFourTensor RankTwoTensor::lkJIDot(const RankTwoTensor &a) const{
    RankFourTensor temp(0.0);
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
                    temp(i,j,k,l)=(*this)(l,k)*a(j,i);
                }
            }
        }
    }
    return temp;
}
RankFourTensor RankTwoTensor::ODot(const RankTwoTensor &a) const{
    // extremely useful for:
    //    dA^-1/dA=-A^-1 \otimes A^-1
    //            =-0.5*Ainv\otimes Ainv=rank-4 tensor
    // for nonlinear constitutive law
    // the proof can be found here:
    // https://en.wikipedia.org/wiki/Tensor_derivative_(continuum_mechanics)

    RankFourTensor temp(0.0);
    for(int i=1;i<=_N;++i){
        for(int j=1;j<=_N;++j){
            for(int k=1;k<=_N;++k){
                for(int l=1;l<=_N;++l){
                    temp(i,j,k,l)=(*this)(i,k)*a(j,l)+(*this)(i,l)*a(j,k);
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
void RankTwoTensor::CalcEigenValueAndEigenVectors(double (&eigval)[3],RankTwoTensor &eigvec)const {
    Eigen::EigenSolver<Eigen::Matrix3d> _eigen_solver;
    Eigen::Matrix3d _M;

    _M<<(*this)(1,1),(*this)(1,2),(*this)(1,3),
        (*this)(2,1),(*this)(2,2),(*this)(2,3),
        (*this)(3,1),(*this)(3,2),(*this)(3,3);
    _eigen_solver.compute(_M);
    eigval[0]=_eigen_solver.eigenvalues()(0).real();
    eigval[1]=_eigen_solver.eigenvalues()(1).real();
    eigval[2]=_eigen_solver.eigenvalues()(2).real();
    for(int i=0;i<_N;++i){
        eigvec(1,i+1)=_eigen_solver.eigenvectors()(0,i).real();
        eigvec(2,i+1)=_eigen_solver.eigenvectors()(1,i).real();
        eigvec(3,i+1)=_eigen_solver.eigenvectors()(2,i).real();
    }
}
//***********************************************
RankFourTensor RankTwoTensor::CalcPostiveProjTensor(double (&eigval)[3],RankTwoTensor &eigvec) const{
    // remember, the eigen vec and eigen value should be used in your material
    // code to calculate the stress and the related constitutive law

    // Algorithm is taken from:
    // C. Miehe and M. Lambrecht, Commun. Numer. Meth. Engng 2001; 17:337~353
    // https://onlinelibrary.wiley.com/doi/epdf/10.1002/cnm.404

    CalcEigenValueAndEigenVectors(eigval,eigvec);

    // C=F^T F=lambda_i M_i where M_i=n_i x n_i
    //         lambda_i-->eigenvalue  n_i--> the eigenvector

    double epos[3],diag[3];
    for(int i=0;i<_N;++i){
        epos[i]=0.5*(abs(eigval[i])+eigval[i]);
        diag[i]=0.0;
        if(eigval[i]>0.0){
            diag[i]=1.0;
        }
    }
    RankFourTensor ProjPos(0.0);
    RankTwoTensor Ma(0.0),Mb(0.0);

    ProjPos.SetToZeros();

    // calculate Ma defined in Eq.(9)-2
    // Ma=n_a x n_a
    for(int i=1;i<=_N;++i){
        Ma.VectorCrossDot(eigvec.IthCol(i),eigvec.IthCol(i));
        // Eq.(19), first term on the right side
        ProjPos+=Ma.CrossDot(Ma)*diag[i-1];
    }

    // Now we calculate the Gab and Gba
    // We need a new rank-2 tensor Mb(same defination as Ma)
    double theta_ab;// defined in Eq.(21)-1
    RankFourTensor Gab(0.0),Gba(0.0);
    const double tol=1.0e-13;
    for(int a=0;a<_N;++a){
        for(int b=0;b<a;++b){
            Ma.VectorCrossDot(eigvec.IthCol(a+1),eigvec.IthCol(a+1));// Eq.(12)
            Mb.VectorCrossDot(eigvec.IthCol(b+1),eigvec.IthCol(b+1));// change the order of Eq.(12)

            Gab=Ma.IkJlDot(Mb)+Ma.IlJkDot(Mb);
            Gba=Mb.IkJlDot(Ma)+Mb.IlJkDot(Ma);
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
