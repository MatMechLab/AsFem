//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.30
//+++ Purpose: Function for reading the whole input file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadInputFile(Mesh &mesh,
                                DofHandler &dofHandler,
                                ElmtSystem &elmtSystem,
                                MateSystem &mateSystem,
                                BCSystem &bcSystem,
                                ICSystem &icSystem,
                                FE &fe){
    ifstream in;
    string str;
    int linenum=0;

    bool HasMeshBlock=false;
    bool HasDofsBlock=false;
    bool HasElmtBlock=false;
    bool HasMateBlock=false;
    bool HasBCBlock=false;
    bool HasICBlock=false;
    bool HasQPointBlock=false;

    if(_HasInputFileName){
        in.open(_InputFileName.c_str(),ios::in);
        while(!in.is_open()){
            MessagePrinter::PrintErrorTxt("can\'t open the input file");
            PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the correct input file name:");
            cin>>_InputFileName;
        }
    }
    else{
        PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the correct input file name:");
            cin>>_InputFileName;
        in.open(_InputFileName.c_str(),ios::in);
        while(!in.is_open()){
            MessagePrinter::PrintErrorTxt("can\'t open the input file");
            PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the input file name:");
            cin>>_InputFileName;
        }
        _HasInputFileName=true;
    }

    linenum=0;

    HasMeshBlock=false;
    HasDofsBlock=false;
    HasElmtBlock=false;
    HasMateBlock=false;
    HasBCBlock=false;
    HasICBlock=false;
    HasQPointBlock=false;

    while(!in.eof()){
        getline(in,str);linenum+=1;
        str=StringUtils::RemoveStrSpace(str);
        str=StringUtils::StrToLower(str);
        if(StringUtils::IsCommentLine(str)||str.size()<1) continue;

        if(str.find("[mesh]")!=string::npos){
            if(!StringUtils::IsBracketMatch(in,linenum)){
                MessagePrinter::PrintErrorTxt("[mesh]/[end] bracket pair dosen\'t match, please check your input file");
                MessagePrinter::AsFem_Exit();
                return false;
            }
            if(ReadMeshBlock(in,str,linenum,mesh)){
                HasMeshBlock=true;
            }
            else{
                HasMeshBlock=false;
            }
        }
        else if(str.find("[dofs]")!=string::npos){
            if(!StringUtils::IsBracketMatch(in,linenum)){
                MessagePrinter::PrintErrorTxt("[dofs]/[end] bracket pair dosen\'t match, please check your input file");
                MessagePrinter::AsFem_Exit();
                return false;
            }
            if(ReadDofsBlock(in,str,linenum,dofHandler)){
                HasDofsBlock=true;
            }
            else{
                HasDofsBlock=false;
            }
        }
        else if(str.find("[elmts]")!=string::npos&&str.find("[ielmts]")==string::npos){
            if(!HasDofsBlock){
                MessagePrinter::PrintErrorTxt("[elmts] block requires the [dofs] block, you should define the [dofs] block first, then given the [elmts] block");
                return false;
            }
            int lastendlinenum;
            if(StringUtils::IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadElmtBlock(in,str,lastendlinenum,linenum,elmtSystem,dofHandler)){
                    HasElmtBlock=true;
                }
                else{
                    HasElmtBlock=false;
                }
            }
            else{
                MessagePrinter::PrintErrorTxt("[elmts]/[end] bracket pair is not match, please check your input file");
                MessagePrinter::AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[mates]")!=string::npos){
            int lastendlinenum;
            if(StringUtils::IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadMateBlock(in,str,lastendlinenum,linenum,mateSystem)){
                    HasMateBlock=true;
                }
                else{
                    HasMateBlock=false;
                }
            }
            else{
                MessagePrinter::PrintErrorTxt("[mates]/[end] bracket pair is not match, please check your input file");
                MessagePrinter::AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[bcs]")!=string::npos){
            if(!HasDofsBlock){
                MessagePrinter::PrintErrorTxt("[bcs] block requires the [dofs] block, you should define the [dofs] block before the [bcs] block");
                MessagePrinter::AsFem_Exit();
                return false;
            }
            int lastendlinenum;
            if(StringUtils::IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadBCBlock(in,str,lastendlinenum,linenum,bcSystem,dofHandler)){
                    HasBCBlock=true;
                }
                else{
                    HasBCBlock=false;
                }
            }
            else{
                MessagePrinter::PrintErrorTxt("[bcs]/[end] bracket pair is not match, please check your input file");
                MessagePrinter::AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[ics]")!=string::npos){
            if(!HasDofsBlock){
                MessagePrinter::PrintErrorTxt("[ics] block requires the [dofs] block, you should define the [dofs] block before the [ics] block");
                MessagePrinter::AsFem_Exit();
                return false;
            }
            int lastendlinenum;
            if(StringUtils::IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadICBlock(in,str,lastendlinenum,linenum,icSystem,dofHandler)){
                    HasICBlock=true;
                }
                else{
                    HasICBlock=false;
                    MessagePrinter::AsFem_Exit();
                }
            }
            else{
                MessagePrinter::PrintErrorTxt("[ics]/[end] bracket pair is not match, please check your input file");
                MessagePrinter::AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[qpoint]")!=string::npos){
            if(!HasMeshBlock){
                MessagePrinter::PrintErrorTxt("[qpoint] block requires the [mesh] block, you should define the [mesh] block first, then given the [qpoint] block");
                return false;
            }
            int lastendlinenum;
            if(StringUtils::IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadQPointBlock(in,str,linenum,fe)){
                    HasQPointBlock=true;
                }
                else{
                    HasQPointBlock=false;
                }
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [qpoint]/[end] bracket pair is not match             !!!   ***\n");
                MessagePrinter::PrintErrorTxt("[qpoint]/[end] bracket pair is not match, please check your input file");
                MessagePrinter::AsFem_Exit();
                return false;
            }
        }
    }

    if(!HasMeshBlock){
        MessagePrinter::PrintErrorTxt("no [mesh] block is found, for FEM analysis, mesh is required");
        MessagePrinter::AsFem_Exit();
    }

    if(!HasDofsBlock){
        MessagePrinter::PrintErrorTxt("no [dofs] block is found, for FEM analysis, dofs are required");
        MessagePrinter::AsFem_Exit();
    }

    if(!HasElmtBlock){
        if(!_IsReadOnly){
            MessagePrinter::PrintErrorTxt("no [elmts] block is found, for FEM analysis, elmts are required");
            MessagePrinter::AsFem_Exit();
        }
    }

    if(!HasMateBlock){
        if(!_IsReadOnly){
            MessagePrinter::PrintWarningTxt("no [mates] block is found, the default material values will be used by the [elmts]",false);
        }
    }

    if(!HasBCBlock){
        if(!_IsReadOnly){
            MessagePrinter::PrintWarningTxt("no [bcs] block is found, no any boundary conditions will be used by AsFem",false);
        }
    }

    if(!HasICBlock){
        if(!_IsReadOnly){
            MessagePrinter::PrintWarningTxt("no [ics] block is found, no any initial conditions will be used by AsFem",false);
        }
    }

    if(!HasQPointBlock){
        if(!_IsReadOnly){
            MessagePrinter::PrintWarningTxt("no [qpoint] block is found, default gauss-legendre integration options will be used by AsFem",false);
        }
        fe.SetDim(mesh.GetDim());
        fe.SetQPointType(QPointType::GAUSSLEGENDRE);
        fe.SetBulkQpOrder(mesh.GetBulkMeshOrder()+1);
        fe.SetBCQpOrder(mesh.GetBulkMeshOrder()+1);
        fe.CreateQPoints(mesh);
    }
    else{
        fe.SetDim(mesh.GetDim());
        fe.CreateQPoints(mesh);
    }

    return true;
}