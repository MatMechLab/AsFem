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
//+++ Date   : 2021.09.03
//+++ Purpose: Calculate the user-defined-material(umat1)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/User1Material.h"

void User1Material::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||Mate.GetScalarMate().size()){}
    Mate.Rank2Materials("stress").SetToZeros();
}

//********************************************************************
void User1Material::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}

    if(InputParams.size()<3){
        MessagePrinter::PrintErrorTxt("for user1material, you need 3 parameters, they are E, nu, and delta");
        MessagePrinter::AsFem_Exit();
    }

    double E0,EE,nu,delta;

    E0=InputParams[1-1];
    nu=InputParams[2-1];
    delta=InputParams[3-1];

    double x,y;
    x=elmtinfo.gpCoords(1);
    y=elmtinfo.gpCoords(2);

    EE=E0*(1.0+delta*sin(x*y));

    RankTwoTensor I,GradU,Strain,Stress,DevStress;
    RankFourTensor Cijkl;

    I.SetToIdentity();

    GradU.SetFromGradU(elmtsoln.gpGradU[1],elmtsoln.gpGradU[2]);
    Strain=(GradU+GradU.Transpose())*0.5;

    Cijkl.SetFromEandNu(EE,nu);

    Stress=Cijkl.DoubleDot(Strain);

    // required by element
    Mate.Rank2Materials("stress")=Stress;
    Mate.Rank2Materials("strain")=Strain;
    Mate.Rank4Materials("jacobian")=Cijkl;

    Mate.ScalarMaterials("MyE")=EE;

    I.SetToIdentity();
    DevStress=Stress-I*(Stress.Trace()/3.0);
    Mate.ScalarMaterials("vonMises")=sqrt(1.5*DevStress.DoubleDot(DevStress));

}




















