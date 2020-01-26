#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::LinearElasticPhaseFieldFracture(const int &isw,const int &nDim,const int &nNodes,const double &JxW,
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
    
    // In your phase field fracture material:
    // _Rank2MateValues[0]--->strain
    // _Rank2MateValues[1]--->stress
    // _Rank2MateValues[2]--->dstress/dD
    // _Rank2MateValues[3]--->dH/dstrain
    // _Rank4MateValues[0]--->elasticity tensor(or your jacobian matrix)
    double Gc=MaterialValues[0];
    double L=MaterialValues[1];
    double viscosity=MaterialValues[2];
    double H=Hist[0];
    double dHdD=Hist[1];
    double valx=0.0,valy=0.0,valz=0.0;
    if(isw==3||isw==6){
        for(int i=1;i<=nNodes;++i){
            if(nDim==2){
                // R_ux
                rhs.coeffRef(3*(i-1)  )+=Rank2MaterialValues[1].IthRow(1)*shp.shape_grad(i)*JxW;
                // R_uy
                rhs.coeffRef(3*(i-1)+1)+=Rank2MaterialValues[1].IthRow(2)*shp.shape_grad(i)*JxW;
                // R_d
                rhs.coeffRef(3*(i-1)+2)+=viscosity*gpV[3-1]*shp.shape_value(i)*JxW
                                        +2*gpU[3-1]*H*shp.shape_value(i)*JxW
                                        -(Gc/(2.0*L))*(1-gpU[3-1])*shp.shape_value(i)*JxW
                                        +2*L*Gc*(gpGradU[3-1]*shp.shape_grad(i))*JxW;

            }
            else if(nDim==3){
                // R_ux
                rhs.coeffRef(4*(i-1)  )+=Rank2MaterialValues[1].IthRow(1)*shp.shape_grad(i)*JxW*ctan[0];
                // R_uy
                rhs.coeffRef(4*(i-1)+1)+=Rank2MaterialValues[1].IthRow(2)*shp.shape_grad(i)*JxW*ctan[0];
                // R_uz
                rhs.coeffRef(4*(i-1)+2)+=Rank2MaterialValues[1].IthRow(3)*shp.shape_grad(i)*JxW*ctan[0];
                // R_c
                rhs.coeffRef(4*(i-1)+3)+=viscosity*gpV[4-1]*shp.shape_value(i)*JxW
                                        +2*gpU[4-1]*H*shp.shape_value(i)*JxW
                                        -(Gc/(2.0*L))*(1-gpU[4-1])*shp.shape_value(i)*JxW
                                        +2*L*Gc*(gpGradU[4-1]*shp.shape_grad(i))*JxW;
            }
            if(isw==6){
                for(int j=1;j<=nNodes;++j){
                    if(nDim==2){
                        // Kux,ux
                        K.coeffRef(3*(i-1)+0,3*(j-1)+0)+=Rank4MaterialValues[0].GetIKjlComponent(1,1,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kux,uy
                        K.coeffRef(3*(i-1)+0,3*(j-1)+1)+=Rank4MaterialValues[0].GetIKjlComponent(1,2,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kux,c
                        K.coeffRef(3*(i-1)+0,3*(j-1)+2)+=shp.shape_value(j)*(Rank2MaterialValues[2].IthRow(1)*shp.shape_grad(i))*JxW*ctan[0];

                        // Kuy,ux
                        K.coeffRef(3*(i-1)+1,3*(j-1)+0)+=Rank4MaterialValues[0].GetIKjlComponent(2,1,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kuy,uy
                        K.coeffRef(3*(i-1)+1,3*(j-1)+1)+=Rank4MaterialValues[0].GetIKjlComponent(2,2,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kuy,c
                        K.coeffRef(3*(i-1)+1,3*(j-1)+2)+=shp.shape_value(j)*(Rank2MaterialValues[2].IthRow(2)*shp.shape_grad(i))*JxW*ctan[0];

                        // Kc,cdot
                        K.coeffRef(3*(i-1)+2,3*(j-1)+2)+=viscosity*shp.shape_value(j)*shp.shape_value(i)*JxW*ctan[1];
                        // Kc,c
                        K.coeffRef(3*(i-1)+2,3*(j-1)+2)+=2*(H+gpU[3-1]*dHdD)*shp.shape_value(j)*shp.shape_value(i)*JxW*ctan[0]
                                                        +(Gc/(2*L))*shp.shape_value(j)*shp.shape_value(i)*JxW*ctan[0]
                                                        +2*Gc*L*(shp.shape_grad(j)*shp.shape_grad(i))*JxW*ctan[0];
                        // Kc,ux
                        valx=0.0;valy=0.0;
                        for(int k=1;k<=3;++k){
                            valx+=0.5*(Rank2MaterialValues[3](1,k)+Rank2MaterialValues[3](k,1))*shp.shape_grad(j).coeff(k-1);
                            valy+=0.5*(Rank2MaterialValues[3](2,k)+Rank2MaterialValues[3](k,2))*shp.shape_grad(j).coeff(k-1);
                        }
                        // Kc,ux
                        K.coeffRef(3*(i-1)+2,3*(j-1)+0)+=2*gpU[3-1]*valx*shp.shape_value(i)*JxW*ctan[0];
                        // Kc,uy
                        K.coeffRef(3*(i-1)+2,3*(j-1)+1)+=2*gpU[3-1]*valy*shp.shape_value(i)*JxW*ctan[0];
                    }
                    else if(nDim==3){
                        // Kux,ux
                        K.coeffRef(4*(i-1)+0,4*(j-1)+0)+=Rank4MaterialValues[0].GetIKjlComponent(1,1,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kux,uy
                        K.coeffRef(4*(i-1)+0,4*(j-1)+1)+=Rank4MaterialValues[0].GetIKjlComponent(1,2,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kux,uz
                        K.coeffRef(4*(i-1)+0,4*(j-1)+2)+=Rank4MaterialValues[0].GetIKjlComponent(1,3,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kux,c
                        K.coeffRef(4*(i-1)+0,4*(j-1)+3)+=shp.shape_value(j)*(Rank2MaterialValues[2].IthRow(1)*shp.shape_grad(i))*JxW*ctan[0];

                        // Kuy,ux
                        K.coeffRef(4*(i-1)+1,4*(j-1)+0)+=Rank4MaterialValues[0].GetIKjlComponent(2,1,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kuy,uy
                        K.coeffRef(4*(i-1)+1,4*(j-1)+1)+=Rank4MaterialValues[0].GetIKjlComponent(2,2,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kuy,uz
                        K.coeffRef(4*(i-1)+1,4*(j-1)+2)+=Rank4MaterialValues[0].GetIKjlComponent(2,3,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kuy,c
                        K.coeffRef(4*(i-1)+1,4*(j-1)+3)+=shp.shape_value(j)*(Rank2MaterialValues[2].IthRow(2)*shp.shape_grad(i))*JxW*ctan[0];

                        // Kuz,ux
                        K.coeffRef(4*(i-1)+2,4*(j-1)+0)+=Rank4MaterialValues[0].GetIKjlComponent(3,1,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kuz,uy
                        K.coeffRef(4*(i-1)+2,4*(j-1)+1)+=Rank4MaterialValues[0].GetIKjlComponent(3,2,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kuz,uy
                        K.coeffRef(4*(i-1)+2,4*(j-1)+2)+=Rank4MaterialValues[0].GetIKjlComponent(3,3,shp.shape_grad(i),shp.shape_grad(j))*JxW*ctan[0];
                        // Kuz,c
                        K.coeffRef(4*(i-1)+2,4*(j-1)+3)+=shp.shape_value(j)*(Rank2MaterialValues[2].IthRow(3)*shp.shape_grad(i))*JxW*ctan[0];

                        // Kc,cdot
                        K.coeffRef(4*(i-1)+3,4*(j-1)+3)+=viscosity*shp.shape_value(j)*shp.shape_value(i)*JxW*ctan[1];
                        // Kc,c
                        K.coeffRef(4*(i-1)+3,4*(j-1)+3)+=2*(H+gpU[4-1]*dHdD)*shp.shape_value(j)*shp.shape_value(i)*JxW*ctan[0]
                                                        +(Gc/(2*L))*shp.shape_value(j)*shp.shape_value(i)*JxW*ctan[0]
                                                        +2*Gc*L*(shp.shape_grad(j)*shp.shape_grad(i))*JxW*ctan[0];
                        // Kc,ux
                        valx=0.0;valy=0.0;valz=0.0;
                        for(int k=1;k<=3;++k){
                            valx+=0.5*(Rank2MaterialValues[3](1,k)+Rank2MaterialValues[3](k,1))*shp.shape_grad(j).coeff(k);
                            valy+=0.5*(Rank2MaterialValues[3](2,k)+Rank2MaterialValues[3](k,2))*shp.shape_grad(j).coeff(k);
                            valz+=0.5*(Rank2MaterialValues[3](3,k)+Rank2MaterialValues[3](k,3))*shp.shape_grad(j).coeff(k);
                        }
                        // Kc,ux
                        K.coeffRef(4*(i-1)+3,4*(j-1)+0)+=2*gpU[4-1]*valx*shp.shape_value(i)*JxW*ctan[0];
                        // Kc,uy
                        K.coeffRef(4*(i-1)+3,4*(j-1)+1)+=2*gpU[4-1]*valy*shp.shape_value(i)*JxW*ctan[0];
                         // Kc,uz
                        K.coeffRef(4*(i-1)+3,4*(j-1)+2)+=2*gpU[4-1]*valz*shp.shape_value(i)*JxW*ctan[0];
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
        //Hist=HistOld;
        if(Hist.size()){};
        if(HistOld.size()) {};
    }
    else if(isw==9){
        // do projection
        Proj[0]=MaterialValues[3];// vonMises stress
        Proj[1]=MaterialValues[4];// hydrostatic stress
        // to get rid of warnings
        Proj[5]=gpGradU[3-1](0);
        Proj[6]=gpGradU[3-1](1);
        Proj[7]=gpU[3-1];
        if(nDim==2){
            Proj[2]=Rank2MaterialValues[1](1,1);//sigxx
            Proj[3]=Rank2MaterialValues[1](2,2);//sigyy
            Proj[4]=Rank2MaterialValues[1](1,2);//sigxy

            //Proj[5]=Rank2MaterialValues[0](1,1);//strxx
            //Proj[6]=Rank2MaterialValues[0](2,2);//stryy
            //Proj[7]=Rank2MaterialValues[0](1,2);//strxy
            Proj[8]=gpU[0];
        }
        else if(nDim==3){
            Proj[2]=Rank2MaterialValues[1](1,1);//sigxx
            Proj[3]=Rank2MaterialValues[1](2,2);//sigyy
            Proj[4]=Rank2MaterialValues[1](3,3);//sigzz
            Proj[5]=Rank2MaterialValues[0](2,3);//sigyz
            Proj[6]=Rank2MaterialValues[0](1,3);//sigxz
            Proj[7]=Rank2MaterialValues[0](1,2);//sigxy
        }
    }
}