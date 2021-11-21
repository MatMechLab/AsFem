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
//+++ Date   : 2021.11.12
//+++ Purpose: Calculate the material properties required by Kobayashi
//+++          element. In this code, we can define:
//+++           1) K
//+++           2) dK
//+++           3) dKdGradEta
//+++           4) ddKdGradEta
//+++ Ref.   :  https://doi.org/10.1016/0167-2789(93)90120-P
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MateSystem/BulkMaterialBase.h"

/**
 * This class calculate the material properties for Kobayashi equation.
 */
class KobayashiMaterial:public BulkMaterialBase{
public:
    /**
     * Initialze the material properties
     */
    virtual void InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate) override;

    /**
     * Calculate the material properties for Kobayashi element
     */
    virtual void ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) override;

private:
    double eta,T;
    double K0,K,dK,ddK,delta,N;
    Vector3d dKdGradEta,ddKdGradEta;
    double dfdeta,d2fdeta2,d2fdetadT;
    double m,dmdT;
    Vector3d GradEta;

    inline double Sign(const double &x){
        return x >= 0.0 ? 1.0 : -1.0;
    }

};
