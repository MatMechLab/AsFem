#include "MaterialSystem/MaterialSystem.h"

void MaterialSystem::LinearElasticPhaseFieldFracMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld){
    
    // just to get rid of warning to error compile mistakes
    // for normal user, you don't need to care about this(warning will always be warning)
    // for developer, please always remember to enable "warning to error" compile flags!!!
    if(nDim){}
    if(t||dt){}
    if(gpCoord(0)){}
    if(gpU[0]){}
    if(gpV[0]){}
    if(gpGradU[0](0)){}
    if(gpGradV[0](0)){}
    if(gpHist[0]){}
    if(gpHistOld[0]){}

    if(InputParams.size()<5){
        cout<<"*** Error: for phase field fracutre,5 params are required!!***"<<endl;
        cout<<"***        you should give E,nu,Gc,L and viscosity !!!     ***"<<endl;
        Msg_AsFem_Exit();
    }
    // InputParams[0]--->E
    // InputParams[1]--->nu
    // InputParams[2]--->Gc
    // InputParams[3]--->L
    // InputParams[4]--->viscosity

    // so in matevalues, we store:
    // Gc, L, viscosity
    // Psi, Psi_pos, Psi_neg
    // Hist,dHistdc(should be 0 currently)
    
    _MateValues[0]=InputParams[2];// for Gc
    _MateValues[1]=InputParams[3];// for L
    _MateValues[2]=InputParams[4];// viscosity

    // here we use Rank4-1 to store the linear Cijkl, the final one should be rank4-0
    _Rank4TensorMateValues[1].SetToZeros();
    _Rank4TensorMateValues[1].SetFromEandNu(InputParams[0],InputParams[1]);
    const double lambda=_Rank4TensorMateValues[1](1,1,2,2);
    const double mu=_Rank4TensorMateValues[1](1,2,1,2);

    RankTwoTensor Stress(0.0),StressPos(0.0),StressNeg(0.0);
    RankTwoTensor Strain(0.0),StrainPos(0.0),StrainNeg(0.0);

    // calculate the elastic strain
    _Rank2TensorMateValues[2].SetToZeros();
    if(nDim==2){
        _Rank2TensorMateValues[2].SetFromGradU(gpGradU[0],gpGradU[1]);
    }
    else{
        _Rank2TensorMateValues[2].SetFromGradU(gpGradU[0],gpGradU[1],gpGradU[2]);
    }
    // this is our strain
    Strain=(_Rank2TensorMateValues[2]+_Rank2TensorMateValues[2].Transpose())*0.5;
    _Rank2TensorMateValues[0]=Strain;
    
    RankTwoTensor eigvec;
    double eigval[3];

    RankFourTensor ProjPos=_Rank2TensorMateValues[0].CalcPostiveProjTensor(eigval,eigvec);
    RankFourTensor I4Sym(RankFourTensor::InitIdentitySymmetric4);
    RankFourTensor ProjNeg=I4Sym-ProjPos;

    StrainPos=ProjPos*Strain;
    StrainNeg=Strain-StrainPos;

    double StrainTrace=Strain.Trace();

    double TrPos= (abs(StrainTrace)+StrainTrace)*0.5;
    double TrNeg=-(abs(StrainTrace)-StrainTrace)*0.5;

    // now we can split the positive and negative stress
    RankTwoTensor I(0.0);
    I.SetToIdentity();// Unity tensor
    StressPos=I*lambda*TrPos+StrainPos*2.0*mu;
    StressNeg=I*lambda*TrNeg+StrainNeg*2.0*mu;

    

    // Finally, we can get the stress as:
    double D=gpU[nDim];// for 2D: ux uy d
                       // for 3D: ux uy uz d
    if(D>1.0-1.0e-3) D=1.0;
    if(D<1.0e-3) D=0.0;
    const double k=1.0e-5;
    _Rank2TensorMateValues[1]=StressPos*(D*D+k)+StressNeg;// Sigma
    _Rank2TensorMateValues[2]=StressPos*2.0*D;        // dSigma/dD

    // now is our constitutive law
    double SignPos,SignNeg;
    SignPos=0.0;
    if(StrainTrace>0.0) SignPos=1.0;
    SignNeg=0.0;
    if(StrainTrace<0.0) SignNeg=1.0;
    _Rank4TensorMateValues[0]=(I.CrossDot(I)*lambda*SignPos+ProjPos*2*mu)*(D*D+k)
                             +(I.CrossDot(I)*lambda*SignNeg+ProjNeg*2*mu);

    // for the fracture free energy
    double Psi,PsiPos,PsiNeg;
    PsiPos=0.5*lambda*TrPos*TrPos+mu*((StrainPos*StrainPos).Trace());
    PsiNeg=0.5*lambda*TrNeg*TrNeg+mu*((StrainNeg*StrainNeg).Trace());
    Psi=D*D*PsiPos+PsiNeg;


    // use the fifth one to calculate the vonMises stress
    _Rank2TensorMateValues[4]=_Rank2TensorMateValues[1];
    // stress deviator tensor
    _Rank2TensorMateValues[4](1,1)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();
    _Rank2TensorMateValues[4](2,2)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();
    _Rank2TensorMateValues[4](3,3)-=(1.0/3.0)*_Rank2TensorMateValues[1].Trace();

    // !!! Remember, the MateValues is already occupied by input params
    _MateValues[0]=InputParams[2];// for Gc
    _MateValues[1]=InputParams[3];// for L
    _MateValues[2]=InputParams[4];// viscosity
    _MateValues[3]=PsiPos;// Positve free energy
    _MateValues[4]=PsiNeg;// negative free energy
    _MateValues[5]=Psi;// free energy

    // Important!!! update the history value
    if(_MateValues[3]>gpHistOld[0]){
        gpHist[0]=_MateValues[3];            //H
        gpHist[1]=_MateValues[3]*2.0*D;      //dH/dD
        _Rank2TensorMateValues[3]=StressPos; // dH/dstrain
    }
    else{
        gpHist[0]=gpHistOld[0];
        gpHist[1]=0.0;
        _Rank2TensorMateValues[3]=0.0;  //dH/dstrain 
    }
    gpHist[1]=0.0; // dH/dD should always equal to 0
    _Rank2TensorMateValues[3]=0.0;  //dH/dstrain 


    // vonMises stress, taken from:
    // https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
    // vonMises=sqrt(1.5*sij*sij)
    // sij is the deviator tensor of stress
    _MateValues[6]=sqrt(1.5*_Rank2TensorMateValues[4].DoubleDot(_Rank2TensorMateValues[4]));
    // hydrostatic stress
    _MateValues[7]=_Rank2TensorMateValues[1].Trace()/3.0;

    if(nDim==2){
        _MateValues[ 8]=_Rank2TensorMateValues[1](1,1);//sigxx
        _MateValues[ 9]=_Rank2TensorMateValues[1](2,2);//sigyy
        _MateValues[10]=_Rank2TensorMateValues[1](1,2);//sigxy

        _MateValues[11]=_Rank2TensorMateValues[0](1,1);//strxx
        _MateValues[12]=_Rank2TensorMateValues[0](2,2);//stryy
        _MateValues[13]=_Rank2TensorMateValues[0](1,2);//strxy
    }
    else if(nDim==3){
        _MateValues[ 8]=_Rank2TensorMateValues[1](1,1);//sigxx
        _MateValues[ 9]=_Rank2TensorMateValues[1](2,2);//sigyy
        _MateValues[10]=_Rank2TensorMateValues[1](3,3);//sigzz
        _MateValues[11]=_Rank2TensorMateValues[1](2,3);//sigyz
        _MateValues[12]=_Rank2TensorMateValues[1](1,3);//sigzx
        _MateValues[13]=_Rank2TensorMateValues[1](1,2);//sigxy

        _MateValues[14]=_Rank2TensorMateValues[0](1,1);//strxx
        _MateValues[15]=_Rank2TensorMateValues[0](2,2);//stryy
        _MateValues[16]=_Rank2TensorMateValues[0](3,3);//strzz
        _MateValues[17]=_Rank2TensorMateValues[0](2,3);//stryz
        _MateValues[18]=_Rank2TensorMateValues[0](1,3);//strzx
        _MateValues[19]=_Rank2TensorMateValues[0](1,2);//strxy
    }
}