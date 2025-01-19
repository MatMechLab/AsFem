//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.11.05
//+++ Purpose: implement the rotated dirichlet boundary condition
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/RotatedDirichletBC.h"

RotatedDirichletBC::RotatedDirichletBC(){
    m_LocalU.resize(11,0.0);// the maximum dofs of each bc block is 10!
}

void RotatedDirichletBC::computeBCValue(const FECalcType &CalcType,
                                        const double &Penalty,
                                        const double &BCValue,
                                        const nlohmann::json &Params,
                                        const LocalElmtInfo &ElmtInfo,
                                        const LocalElmtSolution &ElmtSoln,
                                        const vector<int> &DofIDs,
                                        Vector &U,
                                        SparseMatrix &K,
                                        Vector &RHS){
    if(CalcType==FECalcType::COMPUTERESIDUAL){
        for(int i=0;i<static_cast<int>(DofIDs.size());i++){
            RHS.insertValue(DofIDs[i],0.0);
        }
    }
    else if(CalcType==FECalcType::COMPUTEJACOBIAN){
        for(int i=0;i<static_cast<int>(DofIDs.size());i++){
            K.insertValue(DofIDs[i],DofIDs[i],Penalty);
        }
    }
    else if (CalcType==FECalcType::COMPUTERESIDUALANDJACOBIAN) {
        for(int i=0;i<static_cast<int>(DofIDs.size());i++){
            RHS.insertValue(DofIDs[i],0.0);
            K.insertValue(DofIDs[i],DofIDs[i],Penalty);
        }
    }

    computeU(BCValue,Params,DofIDs,ElmtInfo,ElmtSoln,m_LocalU);
    for(int i=0;i<static_cast<int>(DofIDs.size());i++){
        U.insertValue(DofIDs[i],m_LocalU(i+1));
    }

}

void RotatedDirichletBC::computeU(const double &BCValue,
                                  const nlohmann::json &Params,
                                  const vector<int> &DofIDs,
                                  const LocalElmtInfo &ElmtInfo,
                                  const LocalElmtSolution &ElmtSoln,
                                  VectorXd &LocalU){
    if(ElmtSoln.m_QpU.size()||BCValue){}
    if(ElmtInfo.m_Dim<2){
        MessagePrinter::printErrorTxt("Rotated dirichlet bc works only for 2d surface in a 3d domain");
        MessagePrinter::exitAsFem();
    }
    if(DofIDs.size()<3){
        MessagePrinter::printErrorTxt("Rotated dirichlet must use 3 dofs, namely ux,uz, and uz. Please check your input file");
        MessagePrinter::exitAsFem();
    }
    m_rotation_speed=JsonUtils::getValue(Params,"rotation-speed");
    m_theta=ElmtInfo.m_T*m_rotation_speed*PI/180.0;
    if(JsonUtils::getString(Params,"plane")=="xy"){
        m_x0=JsonUtils::getValue(Params,"x0");
        m_y0=JsonUtils::getValue(Params,"y0");
        m_radius=sqrt((m_x0-ElmtInfo.m_QpCoords0(1))*(m_x0-ElmtInfo.m_QpCoords0(1))
                     +(m_y0-ElmtInfo.m_QpCoords0(2))*(m_y0-ElmtInfo.m_QpCoords0(2)));
        if(m_radius<1.0e-6){
            LocalU(1)=0.0;
            LocalU(2)=0.0;
            LocalU(3)=0.0;
        }
        else{
            if(abs(ElmtInfo.m_QpCoords0(1)-m_x0)<=1.0e-6 && ElmtInfo.m_QpCoords0(2)-m_y0>=0.0){
                // in positive y-axis
                m_theta0=0.5*PI;
            }
            else if(abs(ElmtInfo.m_QpCoords0(1)-m_x0)<=1.0e-6 && ElmtInfo.m_QpCoords0(2)-m_y0<=0.0){
                // in netative y-axis
                m_theta0=-0.5*PI;
            }
            else if(abs(ElmtInfo.m_QpCoords0(2)-m_y0)<=1.0e-6 && ElmtInfo.m_QpCoords0(1)-m_x0>=0.0){
                // in postive x-axis
                m_theta0=0.0;
            }
            else if(abs(ElmtInfo.m_QpCoords0(2)-m_y0)<=1.0e-6 && ElmtInfo.m_QpCoords0(1)-m_x0<=0.0){
                // in negative x-axis
                m_theta0=PI;
            }
            else{
                // for x-y plane
                if(ElmtInfo.m_QpCoords0(1)>m_x0 && ElmtInfo.m_QpCoords0(2)>m_y0){
                    // in 1st phase
                    m_theta0=std::acos((ElmtInfo.m_QpCoords0(1)-m_x0)/m_radius);// in rad
                }
                else if(ElmtInfo.m_QpCoords0(1)<m_x0 && ElmtInfo.m_QpCoords0(2)>m_y0){
                    // in 2nd phase
                    m_theta0=PI-std::acos(-(ElmtInfo.m_QpCoords0(1)-m_x0)/m_radius);
                }
                else if(ElmtInfo.m_QpCoords0(1)<m_x0 && ElmtInfo.m_QpCoords0(2)<m_y0){
                    // in 3rd phase
                    m_theta0=PI+std::acos(-(ElmtInfo.m_QpCoords0(1)-m_x0)/m_radius);
                }
                else if(ElmtInfo.m_QpCoords0(1)>m_x0 && ElmtInfo.m_QpCoords0(2)<m_y0){
                    // in 4th phase
                    m_theta0=2.0*PI-std::acos( (ElmtInfo.m_QpCoords0(1)-m_x0)/m_radius);
                }
            }
            
            LocalU(1)=m_x0+m_radius*std::cos(m_theta0+m_theta)-ElmtInfo.m_QpCoords0(1);
            LocalU(2)=m_y0+m_radius*std::sin(m_theta0+m_theta)-ElmtInfo.m_QpCoords0(2);
            LocalU(3)=0.0;
        }
    }
    else if(JsonUtils::getString(Params,"plane")=="yz"){
        m_y0=JsonUtils::getValue(Params,"y0");
        m_z0=JsonUtils::getValue(Params,"z0");
        m_radius=sqrt((m_y0-ElmtInfo.m_QpCoords0(2))*(m_y0-ElmtInfo.m_QpCoords0(2))
                     +(m_z0-ElmtInfo.m_QpCoords0(3))*(m_z0-ElmtInfo.m_QpCoords0(3)));
        if(m_radius<1.0e-6){
            LocalU(1)=0.0;
            LocalU(2)=0.0;
            LocalU(3)=0.0;
        }
        else{
            if(abs(ElmtInfo.m_QpCoords0(2)-m_y0)<=1.0e-6 && ElmtInfo.m_QpCoords0(3)-m_z0>=0.0){
                // in positive z-axis
                m_theta0=0.5*PI;
            }
            else if(abs(ElmtInfo.m_QpCoords0(2)-m_y0)<=1.0e-6 && ElmtInfo.m_QpCoords0(3)-m_z0<=0.0){
                // in netative z-axis
                m_theta0=-0.5*PI;
            }
            else if(abs(ElmtInfo.m_QpCoords0(3)-m_z0)<=1.0e-6 && ElmtInfo.m_QpCoords0(2)-m_y0>=0.0){
                // in postive y-axis
                m_theta0=0.0;
            }
            else if(abs(ElmtInfo.m_QpCoords0(3)-m_z0)<=1.0e-6 && ElmtInfo.m_QpCoords0(2)-m_y0<=0.0){
                // in negative y-axis
                m_theta0=PI;
            }
            else{
                // for y-z plane
                if(ElmtInfo.m_QpCoords0(2)>m_y0 && ElmtInfo.m_QpCoords0(3)>m_z0){
                    // in 1st phase
                    m_theta0=std::acos((ElmtInfo.m_QpCoords0(2)-m_y0)/m_radius);// in rad
                }
                else if(ElmtInfo.m_QpCoords0(2)<m_y0 && ElmtInfo.m_QpCoords0(3)>m_z0){
                    // in 2nd phase
                    m_theta0=PI-std::acos(-(ElmtInfo.m_QpCoords0(2)-m_y0)/m_radius);
                }
                else if(ElmtInfo.m_QpCoords0(2)<m_y0 && ElmtInfo.m_QpCoords0(3)<m_z0){
                    // in 3rd phase
                    m_theta0=PI+std::acos(-(ElmtInfo.m_QpCoords0(2)-m_y0)/m_radius);
                }
                else if(ElmtInfo.m_QpCoords0(2)>m_y0 && ElmtInfo.m_QpCoords0(3)<m_z0){
                    // in 4th phase
                    m_theta0=2.0*PI-std::acos( (ElmtInfo.m_QpCoords0(2)-m_y0)/m_radius);
                }
            }
            LocalU(1)=0.0;
            LocalU(2)=m_y0+m_radius*std::cos(m_theta0+m_theta)-ElmtInfo.m_QpCoords0(2);
            LocalU(3)=m_z0+m_radius*std::sin(m_theta0+m_theta)-ElmtInfo.m_QpCoords0(3);
        }
    }
    else if(JsonUtils::getString(Params,"plane")=="zx"){
        m_z0=JsonUtils::getValue(Params,"z0");
        m_x0=JsonUtils::getValue(Params,"x0");
        m_radius=sqrt((m_z0-ElmtInfo.m_QpCoords0(3))*(m_z0-ElmtInfo.m_QpCoords0(3))
                     +(m_x0-ElmtInfo.m_QpCoords0(1))*(m_x0-ElmtInfo.m_QpCoords0(1)));
        if(m_radius<1.0e-6){
            LocalU(1)=0.0;
            LocalU(2)=0.0;
            LocalU(3)=0.0;
        }
        else{
            if(abs(ElmtInfo.m_QpCoords0(1)-m_x0)<=1.0e-6 && ElmtInfo.m_QpCoords0(3)-m_z0>=0.0){
                // in positive x-axis
                m_theta0=0.5*PI;
            }
            else if(abs(ElmtInfo.m_QpCoords0(1)-m_x0)<=1.0e-6 && ElmtInfo.m_QpCoords0(3)-m_z0<=0.0){
                // in netative x-axis
                m_theta0=-0.5*PI;
            }
            else if(abs(ElmtInfo.m_QpCoords0(3)-m_z0)<=1.0e-6 && ElmtInfo.m_QpCoords0(1)-m_x0>=0.0){
                // in postive z-axis
                m_theta0=0.0;
            }
            else if(abs(ElmtInfo.m_QpCoords0(3)-m_z0)<=1.0e-6 && ElmtInfo.m_QpCoords0(1)-m_x0<=0.0){
                // in negative z-axis
                m_theta0=PI;
            }
            else{
                // for z-x plane
                if(ElmtInfo.m_QpCoords0(3)>m_z0 && ElmtInfo.m_QpCoords0(1)>m_x0){
                    // in 1st phase
                    m_theta0=std::acos((ElmtInfo.m_QpCoords0(3)-m_z0)/m_radius);// in rad
                }
                else if(ElmtInfo.m_QpCoords0(3)<m_z0 && ElmtInfo.m_QpCoords0(1)>m_x0){
                    // in 2nd phase
                    m_theta0=PI-std::acos(-(ElmtInfo.m_QpCoords0(3)-m_z0)/m_radius);
                }
                else if(ElmtInfo.m_QpCoords0(3)<m_z0 && ElmtInfo.m_QpCoords0(1)<m_x0){
                    // in 3rd phase
                    m_theta0=PI+std::acos(-(ElmtInfo.m_QpCoords0(3)-m_z0)/m_radius);
                }
                else if(ElmtInfo.m_QpCoords0(3)>m_z0 && ElmtInfo.m_QpCoords0(1)<m_x0){
                    // in 4th phase
                    m_theta0=2.0*PI-std::acos( (ElmtInfo.m_QpCoords0(3)-m_z0)/m_radius);
                }
            }
            
            LocalU(3)=m_z0+m_radius*std::cos(m_theta0+m_theta)-ElmtInfo.m_QpCoords0(3);
            LocalU(1)=m_x0+m_radius*std::sin(m_theta0+m_theta)-ElmtInfo.m_QpCoords0(1);
            LocalU(2)=0.0;
        }
    }
    else{
        MessagePrinter::printErrorTxt("The rotated dirichlet bc should be applied to either xy, yz, or zx plane. Please check your input file");
        MessagePrinter::exitAsFem();
    }
}