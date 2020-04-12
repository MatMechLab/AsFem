//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadBCBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,BCSystem &bcSystem,DofHandler &dofHandler){
    // each bc block should looks like:
    //   [bc1]
    //     type=dirichlet [neumann,user1,user2,user3...]
    //     dof=u1
    //     value=1.0 [default is 0.0, so it is not necessary to be given!!!]
    //     boundary=side_name [i.e. left,right]
    //   [end]
    // important: now , str already contains [bcs] !!!

    bool HasBCBlock=false;
    BCBlock bcblock;
    string tempstr,str0;
    vector<double> number;
    vector<string> bclist;

    bool HasElmt=false;
    bool HasValue=false;
    bool HasBoundary=false;
    bool HasDof=false;

    bcblock._IsTimeDependent=false;

    // now str="[bcs]"
    while (linenum<=lastendlinenum){
        getline(in,str);linenum+=1;
        str0=str;
        str=StrToLower(str0);
        str=RemoveStrSpace(str);
        if(IsCommentLine(str)||str.size()<1) continue;

        HasElmt=false;
        HasValue=false;
        HasBoundary=false;
        HasDof=false;

        

        if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
          (str.find("[end]")==string::npos&&str.find("[END]")==string::npos)){
            HasBCBlock=false;
            tempstr=RemoveStrSpace(str);
            if(tempstr.size()<3){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no boundary block name found in [bcs] sub block      !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
            else{
                bcblock.Reset();
                bclist.clear();
                bcblock._BCBlockName=tempstr.substr(1,tempstr.size()-1-1);
                HasBCBlock=true;
            }
            while(str.find("[end]")==string::npos&&str.find("[END]")==string::npos){
                getline(in,str);linenum+=1;
                str0=str;
                str=RemoveStrSpace(str0);
                str=StrToLower(str);
                if(IsCommentLine(str)||str.size()<1) continue;
                if(str.find("type=")!=string::npos){
                    string substr=str.substr(str.find_first_of('=')+1);
                    if(substr.find("dirichlet")!=string::npos && substr.length()==9){
                        bcblock._BCTypeName="dirichlet";
                        bcblock._BCType=BCType::DirichletBC;
                        HasElmt=true;
                    }
                    else if(substr.find("neumann")!=string::npos && substr.length()==7){
                        bcblock._BCTypeName="neumann";
                        bcblock._BCType=BCType::NeumannBC;
                        HasElmt=true;
                    }
                    else if(substr.find("preset")!=string::npos && substr.length()==6){
                        bcblock._BCTypeName="preset";
                        bcblock._BCType=BCType::PresetBC;
                        HasElmt=true;
                    }
                    else if(substr.find("pressure")!=string::npos && substr.length()==8){
                        bcblock._BCTypeName="pressure";
                        bcblock._BCType=BCType::PressureBC;
                        HasElmt=true;
                    }
                    else if(substr.find("user")!=string::npos){
                        number=SplitStrNum(str);
                        if(number.size()<1){
                            Msg_Input_LineError(linenum);
                            PetscPrintf(PETSC_COMM_WORLD,"*** Error: no user number found umat in [bcs] sub block         !!!   ***\n");
                            PetscPrintf(PETSC_COMM_WORLD,"***        type=user1,user2,...,user10 is expected              !!!   ***\n");
                            Msg_AsFem_Exit();
                            return false;
                        }
                        
                        if(int(number[0])<1||int(number[0])>10){
                            Msg_Input_LineError(linenum);
                            PetscPrintf(PETSC_COMM_WORLD,"*** Error: user number is invalid in [bcs] sub block            !!!   ***\n");
                            PetscPrintf(PETSC_COMM_WORLD,"***        type=user1,user2,...,user10 is expected              !!!   ***\n");
                            Msg_AsFem_Exit();
                            return false;
                        }
                        else{
                            bcblock._BCTypeName=str.substr(str.find_first_of('=')+1);
                            switch(int(number[0])){
                                case 1:
                                  bcblock._BCType=BCType::User1BC;
                                  break;
                                case 2:
                                  bcblock._BCType=BCType::User2BC;
                                  break;
                                case 3:
                                  bcblock._BCType=BCType::User3BC;
                                  break;
                                case 4:
                                  bcblock._BCType=BCType::User4BC;
                                  break;
                                case 5:
                                  bcblock._BCType=BCType::User5BC;
                                  break;
                                case 6:
                                  bcblock._BCType=BCType::User6BC;
                                  break;
                                case 7:
                                  bcblock._BCType=BCType::User7BC;
                                  break;
                                case 8:
                                  bcblock._BCType=BCType::User8BC;
                                  break;
                                case 9:
                                  bcblock._BCType=BCType::User9BC;
                                  break;
                                case 10:
                                  bcblock._BCType=BCType::User10BC;
                                  break;
                                default:
                                  bcblock._BCType=BCType::NullBC;
                                  break;
                            }
                            HasElmt=true;
                        }
                    }
                    else{
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported type in [bcs] sub block                  !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        type= %-29s is invalid       !!!   ***\n",substr.c_str());
                        Msg_AsFem_Exit();
                    }
                }
                else if(str.find("dof=")!=string::npos){
                    if(!HasElmt){
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: type= is not given yet in [bcs] sub block            !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        dof= should be given after type= in [bcs] sub block  !!!   ***\n");
                        Msg_AsFem_Exit();
                    }
                    if(str.size()<5){
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: dof name not found in [bcs] sub block                !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        dof=dofname should be given after in [bcs] sub block !!!   ***\n");
                        Msg_AsFem_Exit();
                    }
                    else{
                        bcblock._DofName=str.substr(4);
                        if(!dofHandler.IsValidDofName(bcblock._DofName)){
                            Msg_Input_LineError(linenum);
                            PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid dof name found in [bcs] sub block            !!!   ***\n");
                            PetscPrintf(PETSC_COMM_WORLD,"***        dof must be one of variable in [dofs] block          !!!   ***\n");
                            Msg_AsFem_Exit();
                            return false;
                        }
                        bcblock._DofIndex=dofHandler.GetDofIndexViaName(bcblock._DofName);
                        
                        HasDof=true;
                    }
                }
                else if(str.find("boundary=")!=string::npos){
                    if(!HasElmt){
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: type= is not given yet in [bcs] sub block            !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        boundary= should be given after type= in [bcs] sub block!!!***\n");
                        Msg_AsFem_Exit();
                    }
                    if(str.size()<10){
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: boundary name not found in [bcs] sub block           !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        boundary=bcname should be given after in [bcs] sub block!!!***\n");
                        Msg_AsFem_Exit();
                        return false;
                    }
                    else{
                        int i=str0.find_first_of('=');
                        string substr=str0.substr(i+1,str0.length());
                        bclist=SplitStr(substr,' ');
                        if(bclist.size()>0)
                        {
                            if(IsUniqueStrVec(bclist)){
                                bcblock._BoundaryNameList=bclist;
                                HasBoundary=true;
                            }
                            else{
                                Msg_Input_LineError(linenum);
                                PetscPrintf(PETSC_COMM_WORLD,"*** Error: duplicated boundary name in [bcs] sub block          !!!   ***\n");
                                Msg_AsFem_Exit();
                            }
                        }
                        else{
                            Msg_Input_LineError(linenum);
                            PetscPrintf(PETSC_COMM_WORLD,"*** Error: no boundary name found in [bcs] sub block            !!!   ***\n");
                            PetscPrintf(PETSC_COMM_WORLD,"***        boundary=bc1 bc2 ... is expected                     !!!   ***\n");
                            HasBoundary=false;
                            Msg_AsFem_Exit();
                        }
                    }
                }
                else if(str.find("value=")!=string::npos){
                    if(!HasElmt){
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type= found in [bcs] sub block                    !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        value= must be given after type in [bcs] sub block   !!!   ***\n");
                        Msg_AsFem_Exit();
                    }
                    number=SplitStrNum(str);
                    if(number.size()<1&&str.find("t")==string::npos){
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: no bc value found in [bcs] sub block                 !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        value=real value should be given in [bcs] sub block  !!!   ***\n");
                        Msg_AsFem_Exit();
                        return false;
                    }
                    else{
                        HasValue=true;
                        int i;
                        i=str.find_first_of("=");
                        string substr=str.substr(i+1);
                        // cout<<"bcblock="<<bcblock._BCBlockName<<", str="<<str<<", substr="<<substr<<endl;
                        if(substr.find("t")!=string::npos&&substr.find("*")==string::npos){
                            bcblock._IsTimeDependent=true;
                            bcblock._BCValue=1.0;
                        }
                        else if(substr.find("-t")!=string::npos&&substr.find("*")==string::npos){
                            bcblock._IsTimeDependent=true;
                            bcblock._BCValue=-1.0;
                        }
                        else{
                            if(IsValidExpression(substr)){
                                bcblock._BCValue=number[0];
                                bcblock._IsTimeDependent=true;
                            }
                            else{
                                bcblock._BCValue=number[0];
                                bcblock._IsTimeDependent=false;
                            }
                        }
                    }
                }
            }
            if(HasDof&&HasBoundary&&HasElmt&&HasBCBlock){
                HasBCBlock=true;
                if(!HasValue) bcblock._BCValue=0.0;
                bcSystem.AddBCBlock(bcblock);
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: information is not complete in [mates] sub block     !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        information is missing in [%-20s] !!!   ***\n",bcblock._BCBlockName.c_str());
                if(!HasElmt){
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type found in [bcs] sub block                     !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        type= must be given in each [bcs] sub block          !!!   ***\n");
                }

                if(!HasDof){
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: no dof found in [bcs] sub block                      !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        dof=dof_name must be given in each [bcs] sub block   !!!   ***\n");
                }
               
                if(!HasBoundary){
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: no boundary found in [bcs] sub block                 !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        boundary=bcname must be given in each [bcs] sub block!!!   ***\n");
                }
                return false;
            }
        }
    }

    return true;
}