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
//+++ Date   : 2022.01.29
//+++ Purpose: Calculate the material properties required by thermal
//+++          conduct element. In this code, we can define:
//+++           1) rho--->density
//+++           2) cp --->capacity
//+++           3) qv --->body heat source
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "MateSystem/ThermalMaterial.h"

void ThermalMaterial::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||Mate.GetScalarMate().size()){}

}

//********************************************************************
void ThermalMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}

    if(InputParams.size()<4){
        MessagePrinter::PrintErrorTxt("for thermal material, 4 parameters are required. You need to input: rho, capacity, thermal conductivity, and body heat source");
        MessagePrinter::AsFem_Exit();
    }

    //************************
    //*** here the heat equation is:
    //*** rho*cv*dT/dt=div(k*grad(T))+qv
    //*** 
    
    Mate.ScalarMaterials("rho")=InputParams[0];// density
    Mate.ScalarMaterials("Cp") =InputParams[1];// heat capacity
    Mate.ScalarMaterials("K")  =InputParams[2];// thermal conductivity
    Mate.ScalarMaterials("Q")  =InputParams[3];// body heat source 
    
    // this is designed for general purpose, we keep the derivatives for nonlinear K!
    Mate.ScalarMaterials("dKdT")=0.0;
    Mate.ScalarMaterials("dQdT")=0.0;
    
    Mate.VectorMaterials("gradT")=elmtsoln.gpGradU[1];

}
