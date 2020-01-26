#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::Mechanics(const int &isw,const int &nDim,const int &nNodes,const double &JxW,
                const double &t,const double &dt,const double (&ctan)[2],
                const Eigen::Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &MaterialValues,
                const vector<RankTwoTensor> &Rank2MaterialValues,
                const vector<RankFourTensor> &Rank4MaterialValues,
                vector<double> &Hist,const vector<double> &HistOld,vector<double> &Proj,
                Eigen::MatrixXd &K,Eigen::VectorXd &rhs){

    // just to get rid of complain from compiler(since for developer, warning is treated as errors!)
    // for normal users, you don't need the following two lines,
    // for you, warnings will always be warnings
    // waring will be error only when you are a developer(we suggest developer to do this!!!)
    if(Rank2MaterialValues.size()){};
    if(Rank4MaterialValues.size()){};
    if(gpU[0]){};
    if(gpV[0]){};
    if(gpGradU[0].coeff(0)){};
    if(gpGradV[0].coeff(0)){};
    if(t||dt||gpCoord(0)){};
    
    // In your mechaincal type material:
    // _Rank2MateValues[0]--->strain
    // _Rank2MateValues[1]--->stress
    // _Rank4MateValues[0]--->elasticity tensor(or your jacobian matrix)
    
    if(isw==3||isw==6){
        for(int i=1;i<=nNodes;++i){
            if(nDim==2){
                rhs.coeffRef(2*(i-1)  )+=Rank2MaterialValues[1].IthRow(1)*shp.shape_grad(i)*JxW;
                rhs.coeffRef(2*(i-1)+1)+=Rank2MaterialValues[1].IthRow(2)*shp.shape_grad(i)*JxW;
            }
            else if(nDim==3){
                rhs.coeffRef(3*(i-1)  )+=Rank2MaterialValues[1].IthRow(1)*shp.shape_grad(i)*JxW;
                rhs.coeffRef(3*(i-1)+1)+=Rank2MaterialValues[1].IthRow(2)*shp.shape_grad(i)*JxW;
                rhs.coeffRef(3*(i-1)+2)+=Rank2MaterialValues[1].IthRow(3)*shp.shape_grad(i)*JxW;
            }
            if(isw==6){
                for(int j=1;j<=nNodes;++j){
                    if(nDim==2){
                        // Kux,ux
                        K.coeffRef(2*(i-1)+0,2*(j-1)+0)+=Rank4MaterialValues[0].GetIKjlComponent(1,1,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kux,uy
                        K.coeffRef(2*(i-1)+0,2*(j-1)+1)+=Rank4MaterialValues[0].GetIKjlComponent(1,2,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        
                        // Kuy,ux
                        K.coeffRef(2*(i-1)+1,2*(j-1)+0)+=Rank4MaterialValues[0].GetIKjlComponent(2,1,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kuy,uy
                        K.coeffRef(2*(i-1)+1,2*(j-1)+1)+=Rank4MaterialValues[0].GetIKjlComponent(2,2,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];

                    }
                    else if(nDim==3){
                        // Kux,ux
                        K.coeffRef(3*(i-1)+0,3*(j-1)+0)+=Rank4MaterialValues[0].GetIKjlComponent(1,1,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kux,uy
                        K.coeffRef(3*(i-1)+0,3*(j-1)+1)+=Rank4MaterialValues[0].GetIKjlComponent(1,2,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kux,uz
                        K.coeffRef(3*(i-1)+0,3*(j-1)+2)+=Rank4MaterialValues[0].GetIKjlComponent(1,3,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];

                        // Kuy,ux
                        K.coeffRef(3*(i-1)+1,3*(j-1)+0)+=Rank4MaterialValues[0].GetIKjlComponent(2,1,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kuy,uy
                        K.coeffRef(3*(i-1)+1,3*(j-1)+1)+=Rank4MaterialValues[0].GetIKjlComponent(2,2,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kuy,uy
                        K.coeffRef(3*(i-1)+1,3*(j-1)+2)+=Rank4MaterialValues[0].GetIKjlComponent(2,3,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];

                        // Kuz,ux
                        K.coeffRef(3*(i-1)+2,3*(j-1)+0)+=Rank4MaterialValues[0].GetIKjlComponent(3,1,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kuz,uy
                        K.coeffRef(3*(i-1)+2,3*(j-1)+1)+=Rank4MaterialValues[0].GetIKjlComponent(3,2,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kuz,uy
                        K.coeffRef(3*(i-1)+2,3*(j-1)+2)+=Rank4MaterialValues[0].GetIKjlComponent(3,3,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                    }
                }
            }
        }
    }
    else if(isw==4){
        // init history variable
        fill(Hist.begin(),Hist.end(),0.0);
    }
    else if(isw==8){
        // update hist
        Hist=HistOld;
    }
    else if(isw==9){
        // do projection
        Proj[0]=MaterialValues[0];// vonMises stress
        Proj[1]=MaterialValues[1];// hydrostatic stress
        // to get rid of warnings
        if(nDim==2){
            Proj[2]=Rank2MaterialValues[1](1,1);//sigxx
            Proj[3]=Rank2MaterialValues[1](2,2);//sigxy
            Proj[4]=Rank2MaterialValues[1](1,2);//sigyy

            Proj[5]=Rank2MaterialValues[0](1,1);//strxx
            Proj[6]=Rank2MaterialValues[0](2,2);//strxy
            Proj[7]=Rank2MaterialValues[0](1,2);//stryy
        }
        else if(nDim==3){
            Proj[2]=Rank2MaterialValues[1](1,1);//sigxx
            Proj[3]=Rank2MaterialValues[1](2,2);//sigyy
            Proj[4]=Rank2MaterialValues[1](3,3);//sigzz
            Proj[5]=Rank2MaterialValues[1](2,3);//sigyz
            Proj[6]=Rank2MaterialValues[1](1,3);//sigxz
            Proj[7]=Rank2MaterialValues[1](1,2);//sigxy

            Proj[ 8]=Rank2MaterialValues[0](1,1);//strxx
            Proj[ 9]=Rank2MaterialValues[0](2,2);//stryy
            Proj[10]=Rank2MaterialValues[0](3,3);//strzz
            Proj[11]=Rank2MaterialValues[0](2,3);//stryz
            Proj[12]=Rank2MaterialValues[0](1,3);//strxz
            Proj[13]=Rank2MaterialValues[0](1,2);//strxy
        }
    }
}