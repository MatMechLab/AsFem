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
//+++ Date   : 2022.01.06
//+++ Purpose: Calculate the material properties required by wave 
//+++          propagation equation
//+++          In this code, we can define:
//+++           1) c wave speed
//+++           2) f 
//+++           3) df/du
//+++           4) df/dv 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "MateSystem/WaveMaterial.h"

void WaveMaterial::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||Mate.GetScalarMate().size()){}

}

//********************************************************************
void WaveMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}

    if(InputParams.size()<2){
        MessagePrinter::PrintErrorTxt("for wave material, 3 parameters are required. wave speed, choice of source term, and radius are required");
        MessagePrinter::AsFem_Exit();
    }

    //************************
    //*** here the wave equation is:
    //*** d2u/dt2=c*c*lap(u)+f
    //*** DoFs:
    //***   1-->v(=du/dt)
    //***   2-->u
    
    Mate.ScalarMaterials("C")=InputParams[0];// sigma
    Mate.ScalarMaterials("f")=InputParams[1];// F
    Mate.ScalarMaterials("dfdu")=0.0;// dF/du
    Mate.ScalarMaterials("dfdv")=0.0;

    Mate.VectorMaterials("gradu")=elmtsoln.gpGradU[2];
    int choice=static_cast<int>(InputParams[1]);
    double radius=InputParams[2];
    double x,y,dist;

    x=elmtinfo.gpCoords(1);
    y=elmtinfo.gpCoords(2); 
    dist=sqrt(x*x+y*y); // assume the center is (0,0)
    if(choice==1){
        if(dist<=radius){
            if(elmtinfo.t<0.5){
                Mate.ScalarMaterials("f")=2*sin(x*y);
            }
            else{
                Mate.ScalarMaterials("f")=0.0;
            }
        }
        else{
            Mate.ScalarMaterials("f")=0.0;
        }
    }
    else if(choice==2){
        if(dist<=radius){
            Mate.ScalarMaterials("f")=sin((x*x+y*y)*elmtinfo.t);
        }
        else{
            Mate.ScalarMaterials("f")=0.0;
        }
    }
    else{
        Mate.ScalarMaterials("f")=0.0;
    }

}
