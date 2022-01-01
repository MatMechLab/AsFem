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
//+++ Date   : 2021.02.22
//+++ Purpose: check whether all the settings for each pps is valid
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocess.h"

void Postprocess::CheckWhetherPPSIsValid(const Mesh &mesh){
    for(const auto &pps:_PostProcessBlockList){
        if(pps._PostprocessType==PostprocessType::NODALVALUEPPS){
            if(pps._NodeID<1){
                MessagePrinter::PrintErrorTxt("no nodeid found for Nodevalue postprocess, please check your input file");
                MessagePrinter::AsFem_Exit();
            }
            if(pps._DomainNameList.size()>0||pps._BoundaryNameList.size()>0){
                MessagePrinter::PrintErrorTxt("you dont need domain or side for Nodevalue postprocess, please check your input file");
                MessagePrinter::AsFem_Exit();
            }
        }
        else if(pps._PostprocessType==PostprocessType::AREAPPS){
            for(const auto &side:pps._BoundaryNameList){
                if(mesh.GetBulkMeshDimViaPhyName(side)==mesh.GetBulkMeshDim()){
                    MessagePrinter::PrintErrorTxt("Area postprocess can only be applied for lower dimension sides, not the same dimension as your bulk mesh"
                                                  ", please check your input file");
                    MessagePrinter::AsFem_Exit();
                }
            }
            if(pps._BoundaryNameList.size()<1){
                MessagePrinter::PrintErrorTxt("no boundary or side set name found for your Area postprocess"
                                              ", please check your input file");
                MessagePrinter::AsFem_Exit();
            }
            if(pps._DomainNameList.size()>0){
                MessagePrinter::PrintErrorTxt("you dont need domain for Area postprocess, please check your input file");
                MessagePrinter::AsFem_Exit();
            }
        }
        else if(pps._PostprocessType==PostprocessType::ELEMENTVALUEPPS){
            if(pps._ElementID<1){
                MessagePrinter::PrintErrorTxt("no elmtid found for ElementValue postprocess, please check your input file");
                MessagePrinter::AsFem_Exit();
            }
            if(pps._DomainNameList.size()>0||pps._BoundaryNameList.size()>0){
                MessagePrinter::PrintErrorTxt("you dont need domain or side for ElementValue postprocess, please check your input file");
                MessagePrinter::AsFem_Exit();
            }
        }
        else if(pps._PostprocessType==PostprocessType::VOLUMEPPS){
            for(const auto &domain:pps._DomainNameList){
                if(mesh.GetBulkMeshDimViaPhyName(domain)<mesh.GetBulkMeshDim()){
                    MessagePrinter::PrintErrorTxt("Volume postprocess can only be applied for same dimension domain as your bulk mesh"
                                                  ", please check your input file");
                    MessagePrinter::AsFem_Exit();
                }
            }
            if(pps._DomainNameList.size()<1){
                MessagePrinter::PrintErrorTxt("no domain name found for your Volume postprocess"
                                              ", please check your input file");
                MessagePrinter::AsFem_Exit();
            }
            if(pps._BoundaryNameList.size()>0){
                MessagePrinter::PrintErrorTxt("you dont need side or side for Volume postprocess, please check your input file");
                MessagePrinter::AsFem_Exit();
            }
        }
        else if(pps._PostprocessType==PostprocessType::ELEMENTINTEGRALPPS){
            for(const auto &domain:pps._DomainNameList){
                if(mesh.GetBulkMeshDimViaPhyName(domain)<mesh.GetBulkMeshDim()){
                    MessagePrinter::PrintErrorTxt("ElementalIntegral postprocess can only be applied for same dimension domain as your bulk mesh"
                                                  ", please check your input file");
                    MessagePrinter::AsFem_Exit();
                }
            }
            if(pps._DomainNameList.size()<1){
                MessagePrinter::PrintErrorTxt("no domain name found for your ElementalIntegral postprocess"
                                              ", please check your input file");
                MessagePrinter::AsFem_Exit();
            }
            if(pps._BoundaryNameList.size()>0){
                MessagePrinter::PrintErrorTxt("you dont need side for ElementalIntegral postprocess, please check your input file");
                MessagePrinter::AsFem_Exit();
            }
        }
    }
}