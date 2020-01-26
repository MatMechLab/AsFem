#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::UserElmt1(const int &isw,const int &nDim,const int &nNodes,const double &JxW,
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
    // warning will be error only when you are a developer(we suggest developer to do this!!!)
    if(Rank2MaterialValues.size()){};
    if(Rank4MaterialValues.size()){};
    if(MaterialValues.size()) {};
    if(gpU[0]){};
    if(gpV[0]){};
    if(gpGradU[0].coeff(0)){};
    if(gpGradV[0].coeff(0)){};
    if(t||dt||gpCoord(0)){};

    if(nDim!=2){
        cout<<"*** Error: sorry, uel1 only support for 2D case !!!        ***"<<endl;
        Msg_AsFem_Exit();
    }
    
    // In your umat1,
    // _MateValues[0]--->Young's modulus
    // _MateValues[1]--->Poisson ratio
    // in this uel we don't use any tensor calculation
    // instead, the hand coded one is applied
    const double E0=MaterialValues[0];
    const double nu=MaterialValues[1];

    

    Eigen::Matrix3d C;

    // for plane strain case
    // C(0,0)=1.0-nu;C(0,1)=    nu;C(0,2)=   0.0;
    // C(1,0)=    nu;C(1,1)=1.0-nu;C(1,2)=   0.0;
    // C(2,0)=   0.0;C(2,1)=   0.0;C(2,2)=0.5-nu;
    // C*=E0/((1.0+nu)*(1.0-2.0*nu));

    // for plane stress case
    C(0,0)=1.0;C(0,1)= nu;C(0,2)=         0.0;
    C(1,0)= nu;C(1,1)=1.0;C(1,2)=         0.0;
    C(2,0)=0.0;C(2,1)=0.0;C(2,2)=(1.0-nu)/2.0;
    C*=E0/(1.0-nu*nu);

    // cout<<"C matrix="<<endl<<C<<endl;

    Eigen::Matrix2d F,Ft,E,I;
    Eigen::Vector3d Stress,Strain;

    F(0,0)=1.0+gpGradU[0](0);F(0,1)=    gpGradU[0](1);
    F(1,0)=    gpGradU[1](0);F(1,1)=1.0+gpGradU[1](1);
    // cout<<"F matrix ="<<endl<<F<<endl<<endl;
    // transpose of deformation gradient
    Ft=F.transpose();
    // unity tensor
    I(0,0)=1.0;I(0,1)=0.0;
    I(1,0)=0.0;I(1,1)=1.0;
    E=(Ft*F-I)*0.5;
    

    Strain(0)=E(0,0);Strain(1)=E(1,1);Strain(2)=2.0*E(0,1);
    
    Stress(0)=C(0,0)*Strain(0)+C(0,1)*Strain(1)+C(0,2)*Strain(2);
    Stress(1)=C(1,0)*Strain(0)+C(1,1)*Strain(1)+C(1,2)*Strain(2);
    Stress(2)=C(2,0)*Strain(0)+C(2,1)*Strain(1)+C(2,2)*Strain(2);

    
    if(isw==3||isw==6){
        Eigen::Matrix4d T;
        T.setZero();
        T(0,0)=Stress(0);T(0,1)=Stress(2);
        T(1,0)=Stress(2);T(1,1)=Stress(1);
        T(2,2)=T(0,0);T(2,3)=T(0,1);
        T(3,2)=T(1,0);T(3,3)=T(1,1);
        // for the Linear and nonlinear B matrix
        Eigen::MatrixXd B(3,2*nNodes),Bt(2*nNodes,3);
        Eigen::MatrixXd BNL(4,2*nNodes),BNLt(2*nNodes,4);

        B.setZero();Bt.setZero();
        BNL.setZero();BNLt.setZero();

        for(int i=1;i<=nNodes;i++){
            B(1-1,2*i-1-1)=shp.shape_grad(i)(0)*F(0,0);Bt(2*i-1-1,1-1)=B(1-1,2*i-1-1);
            B(1-1,2*i  -1)=shp.shape_grad(i)(0)*F(1,0);Bt(2*i  -1,1-1)=B(1-1,2*i  -1);

            B(2-1,2*i-1-1)=shp.shape_grad(i)(1)*F(0,1);Bt(2*i-1-1,2-1)=B(2-1,2*i-1-1);
            B(2-1,2*i  -1)=shp.shape_grad(i)(1)*F(1,1);Bt(2*i  -1,2-1)=B(2-1,2*i  -1);

            B(3-1,2*i-1-1)=shp.shape_grad(i)(1)*F(0,0)+shp.shape_grad(i)(0)*F(0,1);
            Bt(2*i-1-1,3-1)=B(3-1,2*i-1-1);
            B(3-1,2*i  -1)=shp.shape_grad(i)(0)*F(1,1)+shp.shape_grad(i)(1)*F(1,0);
            Bt(2*i  -1,3-1)=B(3-1,2*i  -1);

            BNL(1-1,2*i-1-1)=shp.shape_grad(i)(0);BNLt(2*i-1-1,1-1)=BNL(1-1,2*i-1-1);
            BNL(2-1,2*i-1-1)=shp.shape_grad(i)(1);BNLt(2*i-1-1,2-1)=BNL(2-1,2*i-1-1);
            BNL(3-1,2*i  -1)=shp.shape_grad(i)(0);BNLt(2*i  -1,3-1)=BNL(3-1,2*i  -1);
            BNL(4-1,2*i  -1)=shp.shape_grad(i)(1);BNLt(2*i  -1,4-1)=BNL(4-1,2*i  -1);
        }
        Eigen::MatrixXd BtDB(2*nNodes,2*nNodes),BNLtDBNL(2*nNodes,2*nNodes);
        // cout<<"T="<<endl<<T<<endl;

        BtDB=Bt*C*B;
        BNLtDBNL=BNLt*T*BNL;
        // cout<<"BNLtDBNL="<<endl<<BNLtDBNL<<endl;
        for(int i=1;i<=nNodes;++i){
            
            rhs.coeffRef(2*(i-1)  )+=(Bt(2*(i-1)  ,0)*Stress(0)+Bt(2*(i-1)  ,1)*Stress(1)+Bt(2*(i-1)  ,2)*Stress(2))*JxW;
            rhs.coeffRef(2*(i-1)+1)+=(Bt(2*(i-1)+1,0)*Stress(0)+Bt(2*(i-1)+1,1)*Stress(1)+Bt(2*(i-1)+1,2)*Stress(2))*JxW;
            if(isw==6){
                for(int j=1;j<=nNodes;++j){
                    // Kux,ux
                    K.coeffRef(2*(i-1)+0,2*(j-1)+0)+=(BtDB(2*(i-1)+0,2*(j-1)+0)+BNLtDBNL(2*(i-1)+0,2*(j-1)+0))*JxW*ctan[0];
                    // Kux,uy
                    K.coeffRef(2*(i-1)+0,2*(j-1)+1)+=(BtDB(2*(i-1)+0,2*(j-1)+1)+BNLtDBNL(2*(i-1)+0,2*(j-1)+1))*JxW*ctan[0];
                    // Kuy,ux
                    K.coeffRef(2*(i-1)+1,2*(j-1)+0)+=(BtDB(2*(i-1)+1,2*(j-1)+0)+BNLtDBNL(2*(i-1)+1,2*(j-1)+0))*JxW*ctan[0];
                    // Kuy,uy
                    K.coeffRef(2*(i-1)+1,2*(j-1)+1)+=(BtDB(2*(i-1)+1,2*(j-1)+1)+BNLtDBNL(2*(i-1)+1,2*(j-1)+1))*JxW*ctan[0];
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
        Proj[0]=Stress(0);
        Proj[1]=Stress(1);
        Proj[2]=Stress(2);
       
        Proj[3]=Strain(0);//sigxx
        Proj[4]=Strain(1);//sigxy
        Proj[5]=Strain(2);//sigyy
    }
}