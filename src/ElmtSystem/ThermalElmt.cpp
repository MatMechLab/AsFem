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
//+++ Date   : 2021.04.10
//+++ Purpose: implement the residual and jacobian for thermal conduct
//+++          equation.
//+++          rho*cp*dT/dt=div(k*grad(T))+q
//+++          Both k and q could be nonlinear function !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/ThermalElmt.h"

void ThermalElmt::ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
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
        MessagePrinter::PrintErrorTxt("unsupported calculation type in ThermalElmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//***************************************************************************
void ThermalElmt::ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||soln.gpU[0]||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()) {}
    _rho=Mate.ScalarMaterials("rho");
    _K  =Mate.ScalarMaterials("K");
    _Cp =Mate.ScalarMaterials("Cp");
    _Q  =Mate.ScalarMaterials("Q");    
    
    // R_T
    localR(1)=_rho*_Cp*soln.gpV[1]*shp.test
        +_K*(soln.gpGradU[1]*shp.grad_test)
        -_Q*shp.test;

}
//*****************************************************************************
void ThermalElmt::ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU[0]||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()){}

    _rho =Mate.ScalarMaterials("rho");
    _K   =Mate.ScalarMaterials("K");
    _dKdT=Mate.ScalarMaterials("dKdT");
    _Cp  =Mate.ScalarMaterials("Cp");
    _Q   =Mate.ScalarMaterials("Q");
    _dQdT=Mate.ScalarMaterials("dQdT");
    // K_T,T
    localK(1,1)=_rho*_Cp*shp.trial*shp.test*ctan[1]
        +_dKdT*shp.trial*(soln.gpGradU[1]*shp.grad_test)*ctan[0]
        +_K*(shp.grad_trial*shp.grad_test)*ctan[0]
        -_dQdT*shp.trial*shp.test*ctan[0];

}
//*******************************************************************************
void ThermalElmt::ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU[0]||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()||gpProj.size()){}
}
