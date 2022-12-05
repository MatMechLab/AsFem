//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
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
    m_localU.resize(11,0.0);// the maximum dofs of each bc block is 10!
}

void RotatedDirichletBC::computeBCValue(const FECalcType &calctype,const double &penalty,const double &bcvalue,const nlohmann::json &json,
                                const LocalElmtInfo &elmtinfo,
                                const LocalElmtSolution &elmtsoln,
                                const vector<int> &dofids,
                                Vector &U,
                                SparseMatrix &K,
                                Vector &RHS){
    if(calctype==FECalcType::COMPUTERESIDUAL){
        for(int i=0;i<static_cast<int>(dofids.size());i++){
            RHS.insertValue(dofids[i],0.0);
        }
    }
    else if(calctype==FECalcType::COMPUTEJACOBIAN){
        for(int i=0;i<static_cast<int>(dofids.size());i++){
            K.insertValue(dofids[i],dofids[i],penalty);
        }
    }

    computeU(bcvalue,json,dofids,elmtinfo,elmtsoln,m_localU);
    for(int i=0;i<static_cast<int>(dofids.size());i++){
        U.insertValue(dofids[i],m_localU(i+1));
    }

}

void RotatedDirichletBC::computeU(const double &bcvalue,const nlohmann::json &parameters,const vector<int> &dofids,
                          const LocalElmtInfo &elmtinfo,
                          const LocalElmtSolution &elmtsoln,
                          VectorXd &localU){
    if(elmtsoln.m_gpU.size()||bcvalue){}
    if(elmtinfo.m_dim<2){
        MessagePrinter::printErrorTxt("Rotated dirichlet bc works only for 2d surface in a 3d domain");
        MessagePrinter::exitAsFem();
    }
    if(dofids.size()<3){
        MessagePrinter::printErrorTxt("Rotated dirichlet must use 3 dofs, namely ux,uz, and uz. Please check your input file");
        MessagePrinter::exitAsFem();
    }
    m_rotation_speed=JsonUtils::getValue(parameters,"rotation-speed");
    m_theta=elmtinfo.m_t*m_rotation_speed*PI/180.0;
    if(JsonUtils::getString(parameters,"plane")=="xy"){
        m_x0=JsonUtils::getValue(parameters,"x0");
        m_y0=JsonUtils::getValue(parameters,"y0");
        m_radius=sqrt((m_x0-elmtinfo.m_gpCoords0(1))*(m_x0-elmtinfo.m_gpCoords0(1))
                     +(m_y0-elmtinfo.m_gpCoords0(2))*(m_y0-elmtinfo.m_gpCoords0(2)));
        if(m_radius<1.0e-6){
            localU(1)=0.0;
            localU(2)=0.0;
            localU(3)=0.0;
        }
        else{
            if(abs(elmtinfo.m_gpCoords0(1)-m_x0)<=1.0e-6 && elmtinfo.m_gpCoords0(2)-m_y0>=0.0){
                // in positive y-axis
                m_theta0=0.5*PI;
            }
            else if(abs(elmtinfo.m_gpCoords0(1)-m_x0)<=1.0e-6 && elmtinfo.m_gpCoords0(2)-m_y0<=0.0){
                // in netative y-axis
                m_theta0=-0.5*PI;
            }
            else if(abs(elmtinfo.m_gpCoords0(2)-m_y0)<=1.0e-6 && elmtinfo.m_gpCoords0(1)-m_x0>=0.0){
                // in postive x-axis
                m_theta0=0.0;
            }
            else if(abs(elmtinfo.m_gpCoords0(2)-m_y0)<=1.0e-6 && elmtinfo.m_gpCoords0(1)-m_x0<=0.0){
                // in negative x-axis
                m_theta0=PI;
            }
            else{
                // for x-y plane
                if(elmtinfo.m_gpCoords0(1)>m_x0 && elmtinfo.m_gpCoords0(2)>m_y0){
                    // in 1st phase
                    m_theta0=std::acos((elmtinfo.m_gpCoords0(1)-m_x0)/m_radius);// in rad
                }
                else if(elmtinfo.m_gpCoords0(1)<m_x0 && elmtinfo.m_gpCoords0(2)>m_y0){
                    // in 2nd phase
                    m_theta0=PI-std::acos(-(elmtinfo.m_gpCoords0(1)-m_x0)/m_radius);
                }
                else if(elmtinfo.m_gpCoords0(1)<m_x0 && elmtinfo.m_gpCoords0(2)<m_y0){
                    // in 3rd phase
                    m_theta0=PI+std::acos(-(elmtinfo.m_gpCoords0(1)-m_x0)/m_radius);
                }
                else if(elmtinfo.m_gpCoords0(1)>m_x0 && elmtinfo.m_gpCoords0(2)<m_y0){
                    // in 4th phase
                    m_theta0=2.0*PI-std::acos( (elmtinfo.m_gpCoords0(1)-m_x0)/m_radius);
                }
            }
            
            localU(1)=m_x0+m_radius*std::cos(m_theta0+m_theta)-elmtinfo.m_gpCoords0(1);
            localU(2)=m_y0+m_radius*std::sin(m_theta0+m_theta)-elmtinfo.m_gpCoords0(2);
            localU(3)=0.0;
        }
    }
    else if(JsonUtils::getString(parameters,"plane")=="yz"){
        m_y0=JsonUtils::getValue(parameters,"y0");
        m_z0=JsonUtils::getValue(parameters,"z0");
        m_radius=sqrt((m_y0-elmtinfo.m_gpCoords0(2))*(m_y0-elmtinfo.m_gpCoords0(2))
                     +(m_z0-elmtinfo.m_gpCoords0(3))*(m_z0-elmtinfo.m_gpCoords0(3)));
        if(m_radius<1.0e-6){
            localU(1)=0.0;
            localU(2)=0.0;
            localU(3)=0.0;
        }
        else{
            if(abs(elmtinfo.m_gpCoords0(2)-m_y0)<=1.0e-6 && elmtinfo.m_gpCoords0(3)-m_z0>=0.0){
                // in positive z-axis
                m_theta0=0.5*PI;
            }
            else if(abs(elmtinfo.m_gpCoords0(2)-m_y0)<=1.0e-6 && elmtinfo.m_gpCoords0(3)-m_z0<=0.0){
                // in netative z-axis
                m_theta0=-0.5*PI;
            }
            else if(abs(elmtinfo.m_gpCoords0(3)-m_z0)<=1.0e-6 && elmtinfo.m_gpCoords0(2)-m_y0>=0.0){
                // in postive y-axis
                m_theta0=0.0;
            }
            else if(abs(elmtinfo.m_gpCoords0(3)-m_z0)<=1.0e-6 && elmtinfo.m_gpCoords0(2)-m_y0<=0.0){
                // in negative y-axis
                m_theta0=PI;
            }
            else{
                // for y-z plane
                if(elmtinfo.m_gpCoords0(2)>m_y0 && elmtinfo.m_gpCoords0(3)>m_z0){
                    // in 1st phase
                    m_theta0=std::acos((elmtinfo.m_gpCoords0(2)-m_y0)/m_radius);// in rad
                }
                else if(elmtinfo.m_gpCoords0(2)<m_y0 && elmtinfo.m_gpCoords0(3)>m_z0){
                    // in 2nd phase
                    m_theta0=PI-std::acos(-(elmtinfo.m_gpCoords0(2)-m_y0)/m_radius);
                }
                else if(elmtinfo.m_gpCoords0(2)<m_y0 && elmtinfo.m_gpCoords0(3)<m_z0){
                    // in 3rd phase
                    m_theta0=PI+std::acos(-(elmtinfo.m_gpCoords0(2)-m_y0)/m_radius);
                }
                else if(elmtinfo.m_gpCoords0(2)>m_y0 && elmtinfo.m_gpCoords0(3)<m_z0){
                    // in 4th phase
                    m_theta0=2.0*PI-std::acos( (elmtinfo.m_gpCoords0(2)-m_y0)/m_radius);
                }
            }
            localU(1)=0.0;
            localU(2)=m_y0+m_radius*std::cos(m_theta0+m_theta)-elmtinfo.m_gpCoords0(2);
            localU(3)=m_z0+m_radius*std::sin(m_theta0+m_theta)-elmtinfo.m_gpCoords0(3);
        }
    }
    else if(JsonUtils::getString(parameters,"plane")=="zx"){
        m_z0=JsonUtils::getValue(parameters,"z0");
        m_x0=JsonUtils::getValue(parameters,"x0");
        m_radius=sqrt((m_z0-elmtinfo.m_gpCoords0(3))*(m_z0-elmtinfo.m_gpCoords0(3))
                     +(m_x0-elmtinfo.m_gpCoords0(1))*(m_x0-elmtinfo.m_gpCoords0(1)));
        if(m_radius<1.0e-6){
            localU(1)=0.0;
            localU(2)=0.0;
            localU(3)=0.0;
        }
        else{
            if(abs(elmtinfo.m_gpCoords0(1)-m_x0)<=1.0e-6 && elmtinfo.m_gpCoords0(3)-m_z0>=0.0){
                // in positive x-axis
                m_theta0=0.5*PI;
            }
            else if(abs(elmtinfo.m_gpCoords0(1)-m_x0)<=1.0e-6 && elmtinfo.m_gpCoords0(3)-m_z0<=0.0){
                // in netative x-axis
                m_theta0=-0.5*PI;
            }
            else if(abs(elmtinfo.m_gpCoords0(3)-m_z0)<=1.0e-6 && elmtinfo.m_gpCoords0(1)-m_x0>=0.0){
                // in postive z-axis
                m_theta0=0.0;
            }
            else if(abs(elmtinfo.m_gpCoords0(3)-m_z0)<=1.0e-6 && elmtinfo.m_gpCoords0(1)-m_x0<=0.0){
                // in negative z-axis
                m_theta0=PI;
            }
            else{
                // for z-x plane
                if(elmtinfo.m_gpCoords0(3)>m_z0 && elmtinfo.m_gpCoords0(1)>m_x0){
                    // in 1st phase
                    m_theta0=std::acos((elmtinfo.m_gpCoords0(3)-m_z0)/m_radius);// in rad
                }
                else if(elmtinfo.m_gpCoords0(3)<m_z0 && elmtinfo.m_gpCoords0(1)>m_x0){
                    // in 2nd phase
                    m_theta0=PI-std::acos(-(elmtinfo.m_gpCoords0(3)-m_z0)/m_radius);
                }
                else if(elmtinfo.m_gpCoords0(3)<m_z0 && elmtinfo.m_gpCoords0(1)<m_x0){
                    // in 3rd phase
                    m_theta0=PI+std::acos(-(elmtinfo.m_gpCoords0(3)-m_z0)/m_radius);
                }
                else if(elmtinfo.m_gpCoords0(3)>m_z0 && elmtinfo.m_gpCoords0(1)<m_x0){
                    // in 4th phase
                    m_theta0=2.0*PI-std::acos( (elmtinfo.m_gpCoords0(3)-m_z0)/m_radius);
                }
            }
            
            localU(3)=m_z0+m_radius*std::cos(m_theta0+m_theta)-elmtinfo.m_gpCoords0(3);
            localU(1)=m_x0+m_radius*std::sin(m_theta0+m_theta)-elmtinfo.m_gpCoords0(1);
            localU(2)=0.0;
        }
    }
    else{
        MessagePrinter::printErrorTxt("The rotated dirichlet bc should be applied to either xy, yz, or zx plane. Please check your input file");
        MessagePrinter::exitAsFem();
    }
}