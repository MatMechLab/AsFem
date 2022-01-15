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
//+++ Date   : 2022.01.15
//+++ Purpose: implement the cyclic dirichlet boundary condition
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/CyclicDirichletBC.h"

CyclicDirichletBC::CyclicDirichletBC(){
    localU.Resize(11,0.0);// here we set the maximum number of dofs is 10
                      // which means, 'dofs=v1 v2 ... v10', you can append the dofs up to 10!!!
    tspan.resize(50,0.0);// TODO:currently, we assume the time points for the loading is limite to be 50!
    yspan.resize(50,0.0);
} 

void CyclicDirichletBC::ComputeBCValue(const FECalcType &calctype, const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const vector<int> &dofids, const Vector3d &nodecoords, Mat &K, Vec &RHS, Vec &U){
    if(calctype==FECalcType::ComputeResidual){
        for(int i=0;i<static_cast<int>(dofids.size());i++){
            VecSetValue(RHS,dofids[i],0.0,INSERT_VALUES);
        }
    }
    else if(calctype==FECalcType::ComputeJacobian){
        for(int i=0;i<static_cast<int>(dofids.size());i++){
            MatSetValue(K,dofids[i],dofids[i],1.0e16,INSERT_VALUES);
        }
    }
    // the only thing to modify is the 'ComputeU' function!!!
    ComputeU(dofids,bcvalue,params,elmtinfo,nodecoords,localU);
    for(int i=0;i<static_cast<int>(dofids.size());i++){
        VecSetValue(U,dofids[i],localU(i+1),INSERT_VALUES);
    }
}

void CyclicDirichletBC::ComputeU(const vector<int> &dofids, const double &bcvalue, const vector<double> &params, const LocalElmtInfo &elmtinfo, const Vector3d &nodecoords, VectorXd &localU){
    // get rid of unused warnings
    if(bcvalue||params.size()||elmtinfo.dt||nodecoords(1)){}
    // for this simplest case, we just directly apply bcvalue to our localU array
    if(params.size()<4){
        MessagePrinter::PrintErrorTxt("your parameters for cyclic dirichlet bc are too less, at least [x0 y0 x1 y1] should be given");
        MessagePrinter::AsFem_Exit();
    }
    for(int i=0;i<static_cast<int>(params.size()/2);i++){
        tspan[i]=params[2*i];
        yspan[i]=params[2*i+1];
    }
    for(int i=1;i<=static_cast<int>(dofids.size());i++){
        localU(i)=PicewiseLinearInterpolation(tspan,yspan,elmtinfo.t);
    }

}


