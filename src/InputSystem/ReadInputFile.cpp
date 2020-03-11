//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadInputFile(Mesh &mesh,
                                DofHandler &dofHandler,
                                ElmtSystem &elmtSystem,
                                MateSystem &mateSystem,
                                BCSystem &bcSystem,
                                ICSystem &icSystem,
                                Solution &solution,
                                FE &fe,
                                NonlinearSolverBlock &nonlinearSolverBlock,
                                TimeSteppingBlock &timesteppingblock,
                                JobBlock &jobBlock,
                                OutputBlock &outputblock){
    ifstream in;
    string str;
    int linenum=0;

    bool HasMeshBlock=false;
    bool HasDofsBlock=false;
    bool HasBCBlock=false;
    bool HasICBlock=false;
    bool HasElmtBlock=false;
    bool HasMateBlock=false;
    bool HasJobBlock=false;
    bool HasTimeSteppingBlock=false;
    bool HasQpBlock=false;
    // bool HasLinearSolverBlock=false;
    bool HasNonLinearSolverBlock=false;
    bool HasProjectionBlock=false;
    bool HasOutputBlock=false;

    if(_HasInputFileName){
        in.open(_InputFileName.c_str(),ios::in);
        while(!in.is_open()){
            Msg_InputFileName_Invalid(_InputFileName);
            PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the input file name:");
            cin>>_InputFileName;
        }
    }
    else{
        PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the input file name:");
        cin>>_InputFileName;
        in.open(_InputFileName.c_str(),ios::in);
        while(!in.is_open()){
            Msg_InputFileName_Invalid(_InputFileName);
            PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the input file name:");
            cin>>_InputFileName;
        }
        _HasInputFileName=true;
    }

    linenum=0;

    HasMeshBlock=false;
    HasDofsBlock=false;
    HasBCBlock=false;
    HasICBlock=false;
    HasElmtBlock=false;
    HasMateBlock=false;
    HasJobBlock=false;
    HasQpBlock=false;
    // HasLinearSolverBlock=false;
    HasNonLinearSolverBlock=false;
    HasProjectionBlock=false;
    HasTimeSteppingBlock=false;
    HasOutputBlock=false;


    while(!in.eof()){
        getline(in,str);linenum+=1;
        str=RemoveStrSpace(str);
        str=StrToLower(str);
        if(IsCommentLine(str)||str.size()<1) continue;
        if(str.find("[mesh]")!=string::npos){
            if(!IsBracketMatch(in,linenum)){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [mesh]/[end] bracket pair is not match               !!!   ***\n");
                Msg_AsFem_Exit();
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
            if(!IsBracketMatch(in,linenum)){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [dofs]/[end] bracket pair is not match               !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
            if(ReadDofsBlock(in,str,linenum,dofHandler)){
                HasDofsBlock=true;
            }
            else{
                HasDofsBlock=false;
            }
        }
        else if(str.find("[bcs]")!=string::npos){
            if(!HasDofsBlock){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [bcs] block requires [dofs] block                    !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        you should define the [dofs] block before [bcs]      !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadBCBlock(in,str,lastendlinenum,linenum,bcSystem,dofHandler)){
                    HasBCBlock=true;
                }
                else{
                    HasBCBlock=false;
                }
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [bcs]/[end] bracket pair is not match                !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[ics]")!=string::npos){
            if(!HasDofsBlock){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [ics] block requires [dofs] block                    !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        you should define the [dofs] block before [ics]      !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadICBlock(in,str,lastendlinenum,linenum,icSystem,dofHandler)){
                    HasICBlock=true;
                }
                else{
                    HasICBlock=false;
                    Msg_AsFem_Exit();
                }
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [ics]/[end] bracket pair is not match                !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[elmts]")!=string::npos){
            if(!HasDofsBlock){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [elmts] block requires [dofs] block                  !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        you should define the [dofs] block before [elmts]    !!!   ***\n");
                return false;
            }
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadElmtBlock(in,str,lastendlinenum,linenum,elmtSystem,dofHandler)){
                    HasElmtBlock=true;
                }
                else{
                    HasElmtBlock=false;
                }
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [elmts]/[end] bracket pair is not match              !!!   ***\n");
                return false;
            }
        }
        else if(str.find("[mates]")!=string::npos){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadMateBlock(in,str,lastendlinenum,linenum,mateSystem)){
                    HasMateBlock=true;
                }
                else{
                    HasMateBlock=false;
                }
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [mates]/[end] bracket pair is not match              !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[job]")!=string::npos){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadJobBlock(in,str,linenum,jobBlock)){
                    HasJobBlock=true;
                }
                else{
                    HasJobBlock=false;
                }
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [job]/[end] bracket pair is not match                !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[qpoint]")!=string::npos){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadQPointBlock(in,str,linenum,fe)){
                    HasQpBlock=true;
                }
                else{
                    HasQpBlock=false;
                }
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [qpoint]/[end] bracket pair is not match             !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if((str.find("[linearsolver]")!=string::npos)&&str.length()==14){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                // if(ReadLinearSolverBlock(in,str,linenum,linearSolverBlockInfo)){
                //     HasLinearSolverBlock=true;
                // }
                // else{
                //     HasLinearSolverBlock=false;
                // }
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [linearsolver]/[end] bracket pair is not match       !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if((str.find("[nonlinearsolver]")!=string::npos)&&str.length()==17){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadNonlinearSolverBlock(in,str,linenum,nonlinearSolverBlock)){
                    HasNonLinearSolverBlock=true;
                }
                else{
                    HasNonLinearSolverBlock=false;
                }
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [nonlinearsolver]/[end] bracket pair is not match    !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if((str.find("[timestepping]")!=string::npos)&&str.length()==14){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadTimeSteppingBlock(in,str,linenum,timesteppingblock)){
                    HasTimeSteppingBlock=true;
                }
                else{
                    HasTimeSteppingBlock=false;
                }
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [timestepping]/[end] bracket pair is not match       !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if((str.find("[projection]")!=string::npos)&&str.length()==12){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadProjectionBlock(in,str,linenum,solution)){
                    HasProjectionBlock=true;
                }
                else{
                    HasProjectionBlock=false;
                }
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [projection]/[end] bracket pair is not match         !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if((str.find("[output]")!=string::npos)&&str.length()==8){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadOutputBlock(in,str,linenum,outputblock)){
                    HasOutputBlock=true;
                }
                else{
                    HasOutputBlock=false;
                }
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [output]/[end] bracket pair is not match             !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[]")!=string::npos){
            Msg_Input_LineError(linenum);
            Msg_Input_BlockBracketNotComplete();
            Msg_AsFem_Exit();
        }
    }

    if(!HasMeshBlock){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: no [mesh] block found                                !!!   ***\n");
        return false;
    }
    if(!HasDofsBlock){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: no [dofs] block found                                !!!   ***\n");
        return false;
    }

    if(!HasBCBlock){
        PetscPrintf(PETSC_COMM_WORLD,"*** Warning: no [bcs] block found                               !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***          AsFem will treat all boundaries as zero neumann bc !!!   ***\n");
    }

    if(!HasElmtBlock){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: no [elmts] block found                               !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        for a minial input, the [elmts] block is required    !!!   ***\n");
        return false;
    }

    if(!HasMateBlock){
        PetscPrintf(PETSC_COMM_WORLD,"*** Warning: no any materials are given in your input file      !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***          not that case? then please check your input file   !!!   ***\n");
    }

    if(!HasICBlock){
        if(jobBlock._JobType==JobType::TransientJob){
            PetscPrintf(PETSC_COMM_WORLD,"*** Warning: no [ics] block found for a transient analysis      !!!   ***\n");
            PetscPrintf(PETSC_COMM_WORLD,"***          AsFem will treat all dofs as zero                  !!!   ***\n");
        }
    }

    if(!HasQpBlock){
        fe.SetQPointType("gauss");
        fe.SetOrder(mesh.GetMeshOrder()+1);
        fe.SetBCOrder(mesh.GetMeshOrder()+1);
    }


    if(!HasOutputBlock){
        outputblock.Reset();
    }
    outputblock._InputFileName=_InputFileName;

    // if(!HasLinearSolverBlock){
    //     linearSolverBlockInfo.Reset();
    // }

    if(!HasNonLinearSolverBlock){
        nonlinearSolverBlock.Reset();
    }

    if(!HasTimeSteppingBlock){
        timesteppingblock.Reset();
    }

    if(HasProjectionBlock){};// to get rid of compiler's complain

    if(!HasJobBlock){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: no [job] block found                                 !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        for a minial input, the [job] block is required      !!!   ***\n");
        return false;
    }

    //******************************************************************
    //*** now we first set the mate type id for each elmt block
    //******************************************************************
    bool HasMate=false;
    for(PetscInt ielmtblock=1;ielmtblock<=elmtSystem.GetElmtBlocksNum();++ielmtblock){
        if(elmtSystem.GetIthElmtBlock(ielmtblock)._MateBlockName.size()<1){
            elmtSystem.GetIthElmtBlock(ielmtblock)._MateBlockName.clear();
            elmtSystem.GetIthElmtBlock(ielmtblock)._MateType=MateType::NullMate;
            elmtSystem.GetIthElmtBlock(ielmtblock)._MateBlockIndex=1;
        }
        else{
            HasMate=false;
            for(PetscInt imateblock=1;imateblock<=mateSystem.GetMateBlocksNum();++imateblock){
                if(elmtSystem.GetIthElmtBlock(ielmtblock)._MateBlockName==mateSystem.GetIthMateBlock(imateblock)._MateBlockName){
                    HasMate=true;
                    elmtSystem.GetIthElmtBlock(ielmtblock)._MateType=mateSystem.GetIthMateBlock(imateblock)._MateType;
                    elmtSystem.GetIthElmtBlock(ielmtblock)._MateBlockIndex=imateblock;
                }
            }
            if(!HasMate){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no [%-15s] found for [%-15s] in [mates]***\n",elmtSystem.GetIthElmtBlock(ielmtblock)._MateBlockName.c_str(),
                                                                                                                elmtSystem.GetIthElmtBlock(ielmtblock)._ElmtBlockName.c_str());
                PetscPrintf(PETSC_COMM_WORLD,"***        you should define a correct [mates] sub block for it !!!   ***\n");
                Msg_AsFem_Exit();
                return HasMate;
            }
        }
    }
    return true;
}