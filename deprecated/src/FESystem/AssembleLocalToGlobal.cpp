#include "FESystem/FESystem.h"

void FESystem::AssembleLocalToGlobal(const int &isw,const int &ndofs,
                                     const vector<int> &elDofs,
                                     const Eigen::MatrixXd &localK,const Eigen::VectorXd &localR,
                                     Eigen::SparseMatrix<double,Eigen::RowMajor> &AMATRIX,Eigen::VectorXd &RHS){
    if(isw==3||isw==6){
        for(int i=0;i<ndofs;++i){
            RHS.coeffRef(elDofs[i]-1)+=localR.coeffRef(i);
            if(isw==6){
                for(int j=0;j<ndofs;++j){
                    AMATRIX.coeffRef(elDofs[i]-1,elDofs[j]-1)+=-1.0*localK.coeffRef(i,j);
                    if(abs(localK.coeffRef(i,j))>_MaxKMatrixValue) _MaxKMatrixValue=abs(localK.coeffRef(i,j));
                }
            }
        }
    }
}
//************************************
void FESystem::AssembleLocalToGlobal(const int &isw,const int &ndofs,
                                     const vector<int> &elDofs,
                                     const vector<double> &elDofsActiveFlag,
                                     const Eigen::MatrixXd &localK,const Eigen::VectorXd &localR,
                                     Eigen::SparseMatrix<double,Eigen::RowMajor> &AMATRIX,Eigen::VectorXd &RHS){
    if(isw==3||isw==6){
        for(int i=0;i<ndofs;++i){
            // cout<<elDofsActiveFlag[i]<<" ";
            RHS.coeffRef(elDofs[i]-1)+=localR.coeffRef(i)*elDofsActiveFlag[i]*_KMatrixFactor;
            // cout<<localR.coeffRef(i)<<" ";
            if(isw==6){
                for(int j=0;j<ndofs;++j){
                    AMATRIX.coeffRef(elDofs[i]-1,elDofs[j]-1)+=-1.0*localK.coeffRef(i,j)*_KMatrixFactor
                                                               *elDofsActiveFlag[i]*elDofsActiveFlag[j];
                    
                    // _KMatrix factor is used to scaling the K*dx=-R system to have a better convergence behavior
                    // _MaxKMatrixValue is used to apply the dirichlet boundary condition
                    if(abs(localK.coeffRef(i,j))>_MaxKMatrixValue) _MaxKMatrixValue=abs(localK.coeffRef(i,j));
                }
            }
        }
        // cout<<endl;
    }

    // for(int i=0;i<ndofs/2;++i){
    //     for(int j=0;j<ndofs/2;++j){
    //         cout<<localK.coeff(2*i,2*j+1)<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    // for(int i=0;i<ndofs/2;++i){
    //     for(int j=0;j<ndofs/2;++j){
    //         cout<<localK.coeff(2*i+1,2*j)<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<endl;
    
}