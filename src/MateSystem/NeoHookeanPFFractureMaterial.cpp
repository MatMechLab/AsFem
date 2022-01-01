//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.11.07
//+++ Purpose: Calculate the material properties required by Miehe's
//+++          phase field fracture model with neohookean material 
//+++          response
//+++           1) viscosity
//+++           2) Gc
//+++           3) L
//+++           4) H
//+++           5) dHdstrain
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/NeoHookeanPFFractureMaterial.h"

void NeoHookeanPFFractureMaterial::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    //***************************************************
    // get rid of unused warnings
    // **************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||Mate.GetScalarMate().size()){}
    Mate.Rank2Materials("stress").SetToZeros();
    Mate.ScalarMaterials("Hist")=0.0;
}
//***********************************************************
double NeoHookeanPFFractureMaterial::DegradationFun(const double &x) {
    return (1-x)*(1-x);
}
double NeoHookeanPFFractureMaterial::DegradationFunDeriv(const double &x) {
    return 2*(x-1);
}
//************************************************************
void NeoHookeanPFFractureMaterial::ComputeStrain(const LocalElmtInfo &elmtinfo,const LocalElmtSolution &elmtsoln,RankTwoTensor &Strain){

    // here we assume the first dof is d, then the following dofs are ux, uy and uz
    if(elmtinfo.nDim==1){
        MessagePrinter::PrintErrorTxt("Miehe's phase field fracture model only works for 2D and 3D case");
        MessagePrinter::AsFem_Exit();
    }
    else if(elmtinfo.nDim==2){
        // DoFs: d ux uy
        _GradU.SetFromGradU(elmtsoln.gpGradU[2],elmtsoln.gpGradU[3]);
    }
    else if(elmtinfo.nDim==3){
        // DoFs: d ux uy uz
        _GradU.SetFromGradU(elmtsoln.gpGradU[2],elmtsoln.gpGradU[3],elmtsoln.gpGradU[4]);
    }
    _I.SetToIdentity();
    _F=_GradU+_I;
    _E=(_F.Transpose()*_F-_I)*0.5;
    Strain=_E;
}
//************************************************************
void NeoHookeanPFFractureMaterial::ComputeConstitutiveLaws(const vector<double> &InputParams,const double &damage,const RankTwoTensor &strain,RankTwoTensor &Stress,RankFourTensor &Jacobian,const Materials &MateOld, Materials &Mate){
    
    double d;
    double g,dg;// for the degradation function
    const double k=5.0e-5; // for stabilizer

    const double lambda=InputParams[0];
    const double mu=InputParams[1];
    const double Gc=InputParams[2];
    const double L=InputParams[3];
    const double viscosity=InputParams[4];

    Mate.ScalarMaterials("viscosity")=viscosity;
    Mate.ScalarMaterials("Gc")=Gc;
    Mate.ScalarMaterials("L")=L;

    int UseHist=0;
    if(InputParams.size()>=6){
        UseHist=static_cast<int>(InputParams[6-1]);
        if(UseHist<0) UseHist=0;
    }

    if(strain(1,1)){}

    d=damage;
    g= DegradationFun(d);
    dg= DegradationFunDeriv(d);   // derivative of g

    // for our neohookean constitutive laws
    _J=_F.Det();
    _J23=pow(_J,-2.0/3.0);

    _C=_F.Transpose()*_F;
    _Cinv=_C.Inverse();
    
    _I1=_C.Trace();
    _I1bar=_J23*_I1;

    // our elastic free energy density is Psi=U+W
    double U=0.5*lambda*(0.5*(_J*_J-1)-log(_J));
    double W=0.5*mu*(_I1bar-3.0);
    double PsiPos,PsiNeg;

    _StressPos.SetToZeros();
    _StressNeg.SetToZeros();
    _JacPos.SetToZeros();
    _JacNeg.SetToZeros();

    if(_J>=1.0){
        // for tensile loading case
        _I.SetToIdentity();
        PsiPos=U+W;
        PsiNeg=0.0;

        _StressPos=_Cinv*0.5*lambda*(_J*_J-1)
            -_Cinv*(mu/3.0)*_J23*_I1
            +_I*mu*_J23;

        _StressNeg.SetToZeros();

        _JacPos=_Cinv.OTimes(_Cinv)*lambda*_J*_J
            -_Cinv.ODot(_Cinv)*lambda*(_J*_J-1)
            +(
                _Cinv.OTimes(_Cinv)*(_I1/3.0)
               -_Cinv.OTimes(_I)
               +_Cinv.ODot(_Cinv)*_I1
               -_I.OTimes(_Cinv)
             )*(2*mu/3.0)*_J23;

        _JacNeg.SetToZeros();
    }
    else{
        // for compresive loading case
        PsiPos=W;
        PsiNeg=U;
        _I.SetToIdentity();
        _StressPos=_I*mu*_J23-_Cinv*(mu/3.0)*_I1bar;
        _StressNeg=_Cinv*lambda*0.5*(_J*_J-1);

        _JacPos=(
                 _Cinv.OTimes(_Cinv)*(_I1/3.0)
                -_Cinv.OTimes(_I)
                +_Cinv.ODot(_Cinv)*_I1
                -_I.OTimes(_I)
                )*(2*mu/3.0)*_J23;

        _JacNeg=_Cinv.OTimes(_Cinv)*lambda*_J*_J
            -_Cinv.ODot(_Cinv)*lambda*(_J*_J-1);
    }
  
    // store the different elastic free energy densities 
    Mate.ScalarMaterials("Psi")=PsiPos*(g+k)+PsiNeg; 
    Mate.ScalarMaterials("PsiPos")=PsiPos;
    Mate.ScalarMaterials("PsiNeg")=PsiNeg;


    Stress=_StressPos*(g+k)+_StressNeg;
    _PK2=Stress;
    Stress.SetToZeros();
    Stress=_F*_PK2;
    Mate.Rank2Materials("dstressdD")=_F*_StressPos*dg;

    // for jacobian
    Jacobian.SetToZeros();
    I4.SetToIdentity4();
    Jacobian=I4*_PK2+_F*(_JacPos*(g+k)+_JacNeg)*_F.Transpose();

    Mate.ScalarMaterials("H")=0.0;// we use H instead of Hist for our element
                                  // in such a way, users can use the stagger solution!
    if(PsiPos>MateOld.ScalarMaterials("Hist")){
        Mate.ScalarMaterials("Hist")=PsiPos;
        Mate.Rank2Materials("dHdstrain")=_StressPos;
    }
    else{
        Mate.ScalarMaterials("Hist")=MateOld.ScalarMaterials("Hist");
        Mate.Rank2Materials("dHdstrain").SetToZeros();
    }

    Mate.ScalarMaterials("H")=Mate.ScalarMaterials("Hist");

    if(UseHist){
        Mate.ScalarMaterials("H")=MateOld.ScalarMaterials("Hist");
        Mate.Rank2Materials("dHdstrain").SetToZeros();
    }
}
//************************************************************
void NeoHookeanPFFractureMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) {

    //***************************************************************
    // get rid of unused warnings 
    // **************************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU.size()||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}

    if(InputParams.size()<5){
        MessagePrinter::PrintErrorTxt("for Miehe's phase field fracture materials, at least 5 parameters are required, you need to give: lambda, mu, Gc, L, viscosity");
        MessagePrinter::AsFem_Exit();
    }

    // 1st dof: damage
    // 2nd dof: ux
    // 3rd dof: uy
    // 4th dof: uz

    ComputeStrain(elmtinfo,elmtsoln,_Strain);


    ComputeConstitutiveLaws(InputParams,elmtsoln.gpU[1],_Strain,_Stress,_Jacobian,MateOld,Mate);
    // for vonMises stress
    _DevStress=_Stress.Dev();
    // vonMises=sqrt(1.5*sij*sij)
    Mate.ScalarMaterials("vonMises")=sqrt(1.5*_DevStress.DoubleDot(_DevStress));

    Mate.Rank2Materials("strain")=_Strain;
    Mate.Rank2Materials("stress")=_Stress;
    Mate.Rank4Materials("jacobian")=_Jacobian;

    // used in AllenCahn fracture element
    Mate.ScalarMaterials("dFdD")=elmtsoln.gpU[1];
    Mate.ScalarMaterials("d2FdD2")=1.0;
    Mate.ScalarMaterials("M")=1.0/InputParams[5-1];// M=1.0/viscosity

}
