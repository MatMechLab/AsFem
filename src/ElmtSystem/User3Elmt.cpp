//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.09.04
//+++ Purpose: implement the residual and jacobian for the UEL3
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/User3Elmt.h"

void User3Elmt::ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[2],
            const LocalElmtSolution &soln,const LocalShapeFun &shp,
            const Materials &Mate,const Materials &MateOld,
            ScalarMateType &gpProj,
            MatrixXd &localK,VectorXd &localR) {
    if(calctype==FECalcType::ComputeResidual){
        ComputeResidual(elmtinfo,soln,shp,Mate,MateOld,localR);
    }
    else if(calctype==FECalcType::ComputeJacobian){
        ComputeJacobian(elmtinfo,ctan,soln,shp,Mate,MateOld,localK);
    }
    else if(calctype==FECalcType::Projection){
        ComputeProjection(elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj);
    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported calculation type in User3Elmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//***************************************************************************
void User3Elmt::ComputeResidual(const LocalElmtInfo &elmtinfo,
                                const LocalElmtSolution &soln,
                                const LocalShapeFun &shp,
                                const Materials &Mate,const Materials &MateOld,
                                VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||soln.gpU[0]||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()) {}

   localR(1)=0.0;

}
//*****************************************************************************
void User3Elmt::ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[2],
                                const LocalElmtSolution &soln,
                                const LocalShapeFun &shp,
                                const Materials &Mate,const Materials &MateOld,
                                MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU[0]||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()||shp.test){}

    localK(1,1)=0.0;

}
//*******************************************************************************
void User3Elmt::ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[2],
                                  const LocalElmtSolution &soln,
                                  const LocalShapeFun &shp,
                                  const Materials &Mate,const Materials &MateOld,
                                  ScalarMateType &gpProj) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU[0]||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()||gpProj.size()){}

}
