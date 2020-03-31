//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadICBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,ICSystem &icSystem,DofHandler &dofHandler){
    // each bc block should looks like:
    //   [ic1]
    //     type=const [random,user1,user2,user3...]
    //     dof=u1
    //     params=1.0 [default is 0.0, so it is not necessary to be given!!!]
    //     block=block_name [i.e. all]
    //   [end]
    // important: now , str already contains [bcs] !!!

    bool HasICBlock=false;
    ICBlock icblock;
    string tempstr,str0,stro;
    vector<double> number;

    bool HasElmt=false;
    bool HasValue=false;
    bool HasBlock=false;
    bool HasDof=false;

    // now str="[bcs]"
    while (linenum<=lastendlinenum){
        getline(in,str);linenum+=1;
        str0=str;
        str=RemoveStrSpace(str);
        stro=str;
        str=StrToLower(str);
        if(IsCommentLine(str)||str.size()<1) continue;

        HasElmt=false;
        HasValue=false;
        HasBlock=false;
        HasDof=false;

        icblock.Reset();

        if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
          (str.find("[end]")==string::npos&&str.find("[END]")==string::npos)){
            HasICBlock=false;
            str=StrToLower(str);
            tempstr=RemoveStrSpace(str);
            if(tempstr.size()<3){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no IC block name found in [ics] sub block            !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
            else{
                icblock.Reset();
                icblock._ICBlockName=tempstr.substr(1,tempstr.size()-1-1);
                HasICBlock=true;
            }
            while(str.find("[end]")==string::npos&&str.find("[END]")==string::npos){
                getline(in,str);linenum+=1;
                str0=str;
                str=RemoveStrSpace(str);
                stro=str;
                str=StrToLower(stro);
                if(IsCommentLine(str)||str.size()<1) continue;
                if(str.find("type=")!=string::npos){
                    string substr=str.substr(str.find('=')+1);
                    if(substr.find("const")!=string::npos&&substr.length()==5){
                        icblock._ICTypeName="const";
                        icblock._ICType=ICType::ConstIC;
                        HasElmt=true;
                    }
                    else if(substr.find("random")!=string::npos&&substr.length()==6){
                        icblock._ICTypeName="random";
                        icblock._ICType=ICType::RandomIC;
                        HasElmt=true;
                    }
                    else if(substr.find("circle")!=string::npos&&substr.length()==6){
                        icblock._ICTypeName="circle";
                        icblock._ICType=ICType::CircleIC;
                        HasElmt=true;
                    }
                    else if(substr.find("rectangle")!=string::npos&&substr.length()==9){
                        icblock._ICTypeName="rectangle";
                        icblock._ICType=ICType::RectangleIC;
                        HasElmt=true;
                    }
                    else if(substr.find("user")!=string::npos){
                        number=SplitStrNum(str);
                        if(number.size()<1){
                            Msg_Input_LineError(linenum);
                            PetscPrintf(PETSC_COMM_WORLD,"*** Error: no user number found umat in [ics] sub block         !!!   ***\n");
                            PetscPrintf(PETSC_COMM_WORLD,"***        type=user1,user2,...,user10 is expected              !!!   ***\n");
                            Msg_AsFem_Exit();
                            return false;
                        }
                        
                        if(int(number[0])<1||int(number[0])>10){
                            Msg_Input_LineError(linenum);
                            PetscPrintf(PETSC_COMM_WORLD,"*** Error: user number is invalid in [ics] sub block            !!!   ***\n");
                            PetscPrintf(PETSC_COMM_WORLD,"***        type=user1,user2,...,user10 is expected              !!!   ***\n");
                            Msg_AsFem_Exit();
                            return false;
                        }
                        else{
                            icblock._ICTypeName=str.substr(str.find_first_of('=')+1);
                            switch(int(number[0])){
                                case 1:
                                    icblock._ICType=ICType::User1IC;
                                    break;
                                case 2:
                                    icblock._ICType=ICType::User2IC;
                                    break;
                                case 3:
                                    icblock._ICType=ICType::User3IC;
                                    break;
                                case 4:
                                    icblock._ICType=ICType::User4IC;
                                    break;
                                case 5:
                                    icblock._ICType=ICType::User5IC;
                                    break;
                                case 6:
                                    icblock._ICType=ICType::User6IC;
                                    break;
                                case 7:
                                    icblock._ICType=ICType::User7IC;
                                    break;
                                case 8:
                                    icblock._ICType=ICType::User8IC;
                                    break;
                                case 9:
                                    icblock._ICType=ICType::User9IC;
                                    break;
                                case 10:
                                    icblock._ICType=ICType::User10IC;
                                    break;
                                default:
                                    icblock._ICType=ICType::NullIC;
                                    break;
                            }
                            HasElmt=true;
                        }
                    }
                    else{
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported IC type in [ics] sub block               !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        type= %-25s is invalid           !!!   ***\n",substr.c_str());
                        Msg_AsFem_Exit();
                    }
                }
                else if(str.find("dof=")!=string::npos){
                    if(!HasElmt){
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: type= is not given yet in [ics] sub block            !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        dof= should be given after type= in [ics] sub block  !!!   ***\n");
                        Msg_AsFem_Exit();
                    }
                    if(str.size()<5){
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: no dof name found in [ics] sub block                 !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        dof=dof_name should be given in [ics] sub block      !!!   ***\n");
                        Msg_AsFem_Exit();
                    }
                    else{
                        icblock._DofName=stro.substr(4);
                        if(!dofHandler.IsValidDofName(icblock._DofName)){
                            Msg_Input_LineError(linenum);
                            PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid dof name found in [ics] sub block            !!!   ***\n");
                            PetscPrintf(PETSC_COMM_WORLD,"***        dof must be one of variable in [dofs] block          !!!   ***\n");
                            Msg_AsFem_Exit();
                            return false;
                        }
                        icblock._DofIndex=dofHandler.GetDofIndexViaName(icblock._DofName);
                        HasDof=true;
                    }
                }
                else if(str.find("domain=")!=string::npos){
                    if(!HasElmt){
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: type= is not given yet in [ics] sub block            !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        domoain= should be given after type= in [ics] sub block !!!***\n");
                        Msg_AsFem_Exit();
                    }
                    if(str.size()<7){
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: no domain name found in [ics] sub block              !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        domain=domain_name should be given in [ics] sub block!!!   ***\n");
                        Msg_AsFem_Exit();
                        return false;
                    }
                    else{
                        icblock._DomainName=str.substr(str.find_first_of("=")+1);
                        HasBlock=true;
                    }
                }
                else if(str.find("params=")!=string::npos){
                    if(!HasElmt){
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: type= is not given yet in [ics] sub block            !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        params= should be given after type in [ics] sub block!!!   ***\n");
                        Msg_AsFem_Exit();
                    }
                    number=SplitStrNum(str0);
                    if(number.size()<1){
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: no bc params found in [bcs] sub block                !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        params=real value should be given in [bcs] sub block !!!   ***\n");
                        Msg_AsFem_Exit();
                        return false;
                    }
                    else{
                        icblock._Params=number;
                        HasValue=true;
                    }
                }
            }
            if(HasDof&&HasElmt&&HasICBlock){
                HasICBlock=true;
                if(!HasValue) {
                    icblock._Params.clear();
                    icblock._Params.push_back(0.0);
                    icblock._Params.push_back(1.0);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Warning: no params are given in [ics] sub block             !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***          default value will be used in your IC element      !!!   ***\n");
                }
                if(!HasBlock){
                    icblock._DomainName="alldomain";
                    HasBlock=true;
                }
                icSystem.AddICBlock(icblock);
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: information is not complete in [ics] sub block       !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        information is missing in [%-20s] !!!   ***\n",icblock._ICBlockName.c_str());
                if(!HasElmt){
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type found in [ics] sub block                     !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        type= must be given in each [ics] sub block          !!!   ***\n");
                }

                if(!HasDof){
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: no dof found in [ics] sub block                      !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        dof=dof_name must be given in each [ics] sub block   !!!   ***\n");
                }
                return false;
            }
        }
    }

    return true;
}