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
//+++ Date   : 2021.02.21
//+++ Purpose: check if the applied boundary name is consistent with
//+++          the one we defined in our mesh
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

bool BCSystem::CheckAppliedBCNameIsValid(const Mesh &mesh){
    string bctypename;
    vector<string> bcnamelist;
    BCType bctype;
    bool HasFoundName;
    for(const auto &bcblock:_BCBlockList){
        bcnamelist=bcblock._BoundaryNameList;
        bctype=bcblock._BCType;
        HasFoundName=false;
        if(bctype==BCType::NEUMANNBC||
           bctype==BCType::DIRICHLETBC||
           bctype==BCType::PRESSUREBC){
            for(const auto &bcname:bcnamelist){
                HasFoundName=false;
                for(int i=1;i<=mesh.GetBulkMeshPhysicalGroupNum();i++){
                    if(mesh.GetBulkMeshIthPhysicalName(i)==bcname){
                        HasFoundName=true;
                        break;
                    }
                }
                if(!HasFoundName){
                    MessagePrinter::PrintErrorTxt("you applied bc on '"+bcname+"', however, your mesh dosent have this boundary element set "
                                                      "you should apply your elemental type boundary condition to elemental set, instead of node set");
                    return false;
                }
            }
        }
        else if(bctype==BCType::NODALDIRICHLETBC||
                bctype==BCType::NODALNEUMANNBC||
                bctype==BCType::NODALFLUXBC||
                bctype==BCType::NODALFORCEBC){
            for(const auto &bcname:bcnamelist){
                HasFoundName=false;
                for(int i=1;i<=mesh.GetBulkMeshNodeSetPhysicalGroupNum();i++){
                    if(mesh.GetBulkMeshIthNodeSetPhysicalName(i)==bcname){
                        HasFoundName=true;
                        break;
                    }
                }
                if(!HasFoundName){
                    MessagePrinter::PrintErrorTxt("you applied bc on '"+bcname+"', however, your mesh dosent have this boundary node set "
                                                      "you should apply your nodal type boundary condition (type=nodalxxx) to node set, instead of element set");
                    return false;
                }
            }
        }
        else{
            MessagePrinter::PrintWarningTxt("we can not check the user defined bc, "
                                            "you should make sure by yourselves that the boundary name "
                                            "is correct(either using node set or element set)");
        }
    }
    return true;
}