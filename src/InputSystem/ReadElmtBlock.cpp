//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadElmtBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,ElmtSystem &elmtSystem,DofHandler &dofHandler){
    // elmt block format
    // [elmt]
    //  [solid]
    //    type=solid2d
    //    dofs=u1 u2
    //    mate=mate1 [can be ignored]
    //    block=all  [can be ignored]
    //  [end]
    // [end]
    bool HasElmtBlock=false;
    bool HasElmtType=false;
    bool HasDofs=false;
    bool HasBlock=false;
    bool HasMate=false;
    ElmtBlock elmtBlock;
    string tempstr,str0,substr;
    vector<double> number;
    vector<string> strlist;


    elmtBlock.Reset();

    // now the str is [elmts]
    while(linenum<=lastendlinenum){
        getline(in,str);linenum+=1;
        str0=str;
        str=StrToLower(str0);
        str=RemoveStrSpace(str);
        if(IsCommentLine(str)||str.size()<1) continue;

        HasElmtBlock=false;
        HasElmtType=false;
        HasDofs=false;
        HasBlock=false;
        HasMate=false;

        if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
          (str.find("[end]")==string::npos&&str.find("[END]")==string::npos)){
            HasElmtBlock=false;
            tempstr=RemoveStrSpace(str);
            if(tempstr.size()<3){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no elmt block name found in [elmts]                  !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        a sub block like [xx] should be given inside [elmts] !!!   ***\n");
                Msg_AsFem_Exit();
                return false;
            }
            else{
                elmtBlock._ElmtBlockName=tempstr.substr(tempstr.find_first_of('[')+1,tempstr.find_first_of(']')-1);
                HasElmtBlock=true;
            }
            while(str.find("[end]")==string::npos&&str.find("[END]")==string::npos){
                getline(in,str);linenum+=1;
                str0=str;
                str=StrToLower(str0);
                str=RemoveStrSpace(str);

                if(IsCommentLine(str)||str.size()<1) {
                    //cout<<"str="<<str<<endl;
                    continue;
                }

                if(str.find("type=")!=string::npos){
                    substr=str.substr(str.find_first_of('=')+1);
                    if(substr.find("poisson")!=string::npos && substr.length()==7){
                        elmtBlock._ElmtTypeName="poisson";
                        elmtBlock._ElmtType=ElmtType::PoissonElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("diffusion")!=string::npos && substr.length()==9){
                        elmtBlock._ElmtTypeName="diffusion";
                        elmtBlock._ElmtType=ElmtType::DiffusionElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("wave")!=string::npos && substr.length()==4){
                        elmtBlock._ElmtTypeName="wave";
                        elmtBlock._ElmtType=ElmtType::WaveElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("mechanics")!=string::npos && substr.length()==9){
                        elmtBlock._ElmtTypeName="mechanics";
                        elmtBlock._ElmtType=ElmtType::MechanicsElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("cahnhilliard")!=string::npos && substr.length()==12){
                        elmtBlock._ElmtTypeName="cahnhilliard";
                        elmtBlock._ElmtType=ElmtType::CahnHilliardElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("tensorcahnhilliard")!=string::npos && substr.length()==18){
                        elmtBlock._ElmtTypeName="tensorcahnhilliard";
                        elmtBlock._ElmtType=ElmtType::TensorCahnHilliardElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("mechcahnhilliard")!=string::npos && substr.length()==16){
                        elmtBlock._ElmtTypeName="mechcahnhilliard";
                        elmtBlock._ElmtType=ElmtType::MechCahnHilliardElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("thermalmechanics")!=string::npos && substr.length()==16){
                        elmtBlock._ElmtTypeName="thermalmechanics";
                        elmtBlock._ElmtType=ElmtType::ThermalMechanicsElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("phasefieldfracture")!=string::npos && substr.length()==18){
                        elmtBlock._ElmtTypeName="phasefieldfracture";
                        elmtBlock._ElmtType=ElmtType::PhaseFieldFracElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("miehefrac")!=string::npos && substr.length()==9){
                        elmtBlock._ElmtTypeName="miehefrac";
                        elmtBlock._ElmtType=ElmtType::MieheFractureElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("cohesivepffrac")!=string::npos && substr.length()==14){
                        elmtBlock._ElmtTypeName="cohesivepffrac";
                        elmtBlock._ElmtType=ElmtType::CohesivePFFracElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("bordenfrac")!=string::npos && substr.length()==10){
                        elmtBlock._ElmtTypeName="bordenfrac";
                        elmtBlock._ElmtType=ElmtType::BordenFractureElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("thermal")!=string::npos && substr.length()==7){
                        elmtBlock._ElmtTypeName="thermal";
                        elmtBlock._ElmtType=ElmtType::ThermalElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("dendrite")!=string::npos && substr.length()==8){
                        elmtBlock._ElmtTypeName="dendrite";
                        elmtBlock._ElmtType=ElmtType::DendriteElmt;
                        HasElmtType=true;
                    }
                    else if(substr.find("user")!=string::npos){
                        number=SplitStrNum(substr);
                        if(number.size()<1){
                            Msg_Input_LineError(linenum);
                            PetscPrintf(PETSC_COMM_WORLD,"*** Error: no user element number found in [elmts] block        !!!   ***\n");
                            PetscPrintf(PETSC_COMM_WORLD,"***        type=user1,user2,user3... is expected                !!!   ***\n");
                            Msg_AsFem_Exit();
                            return false;
                        }
                        if(int(number[0])<1||int(number[0])>20){
                            Msg_Input_LineError(linenum);
                            PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid user element number in [elmts] block         !!!   ***\n");
                            PetscPrintf(PETSC_COMM_WORLD,"***        type=user1~user20  is expected                       !!!   ***\n");
                            Msg_AsFem_Exit();
                            return false;
                        }
                        else{
                            int iuel=int(number[0]);
                            elmtBlock._ElmtTypeName=str.substr(str.find_first_of('=')+1);
                            if(iuel==1){
                                elmtBlock._ElmtType=ElmtType::User1Elmt;
                            }
                            else if(iuel==2){
                                elmtBlock._ElmtType=ElmtType::User2Elmt;
                            }
                            else if(iuel==3){
                                elmtBlock._ElmtType=ElmtType::User3Elmt;
                            }
                            else if(iuel==4){
                                elmtBlock._ElmtType=ElmtType::User4Elmt;
                            }
                            else if(iuel==5){
                                elmtBlock._ElmtType=ElmtType::User5Elmt;
                            }
                            else if(iuel==6){
                                elmtBlock._ElmtType=ElmtType::User6Elmt;
                            }
                            else if(iuel==7){
                                elmtBlock._ElmtType=ElmtType::User7Elmt;
                            }
                            else if(iuel==8){
                                elmtBlock._ElmtType=ElmtType::User8Elmt;
                            }
                            else if(iuel==9){
                                elmtBlock._ElmtType=ElmtType::User9Elmt;
                            }
                            else if(iuel==10){
                                elmtBlock._ElmtType=ElmtType::User10Elmt;
                            }
                            else if(iuel==11){
                                elmtBlock._ElmtType=ElmtType::User11Elmt;
                            }
                            else if(iuel==12){
                                elmtBlock._ElmtType=ElmtType::User12Elmt;
                            }
                            else if(iuel==13){
                                elmtBlock._ElmtType=ElmtType::User13Elmt;
                            }
                            else if(iuel==14){
                                elmtBlock._ElmtType=ElmtType::User14Elmt;
                            }
                            else if(iuel==15){
                                elmtBlock._ElmtType=ElmtType::User15Elmt;
                            }
                            else if(iuel==16){
                                elmtBlock._ElmtType=ElmtType::User16Elmt;
                            }
                            else if(iuel==17){
                                elmtBlock._ElmtType=ElmtType::User17Elmt;
                            }
                            else if(iuel==18){
                                elmtBlock._ElmtType=ElmtType::User18Elmt;
                            }
                            else if(iuel==19){
                                elmtBlock._ElmtType=ElmtType::User19Elmt;
                            }
                            else if(iuel==20){
                                elmtBlock._ElmtType=ElmtType::User20Elmt;
                            }
                            HasElmtType=true;
                        }
                    }
                    else{
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported element type in [elmts] block            !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        type=%35s is invalid                !!!   ***\n",substr.c_str());
                        Msg_AsFem_Exit();
                        return false;
                    }
                }
                else if(str.find("dofs=")!=string::npos){
                    if(str.size()<5){
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: no dof name found in [elmts] block                   !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        dofs=dof1 dof2 ... is expected                       !!!   ***\n");
                        Msg_AsFem_Exit();
                    }
                    else
                    {
                        substr=str0.substr(str0.find_first_of('=')+1);
                        vector<string> strtmp=SplitStr(substr,' ');
                        if(!dofHandler.IsValidDofNameVec(strtmp)){
                            PetscPrintf(PETSC_COMM_WORLD,"*** Errors: dof name is not valid in [elmts] sub block          !!!   ***\n");
                            PetscPrintf(PETSC_COMM_WORLD,"***         the name should be one list in [dofs] block         !!!   ***\n");
                            return false;
                        }
                        if(IsUniqueStrVec(strtmp)){
                            elmtBlock._DofNameList=SplitStr(substr,' ');
                            elmtBlock._DofIndexList=dofHandler.GetDofsIndexFromNameVec(elmtBlock._DofNameList);
                            HasDofs=true;
                        }
                        else{
                            HasDofs=false;
                            Msg_Input_LineError(linenum);
                            PetscPrintf(PETSC_COMM_WORLD,"*** Error: duplicate dof name found in [elmts] block            !!!   ***\n");
                            PetscPrintf(PETSC_COMM_WORLD,"***        dofs=dof1 dof2 ... is expected                       !!!   ***\n");
                            Msg_AsFem_Exit();
                        }
                    }
                }
                else if(str.find("mate=")!=string::npos){
                    if(str.size()<5){
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: material block name not found in [elmts] block       !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        mate=name_of_material_block is expected              !!!   ***\n");
                        Msg_AsFem_Exit();
                        return false;
                    }
                    else{
                        elmtBlock._MateBlockName=str.substr(str.find_first_of("=")+1);
                        HasMate=true;
                    }
                }
                else if(str.find("domain=")!=string::npos){
                    if(str.size()<7){
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: no domain name found in [elmts] block                !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        domain=name_of_domain is expected                    !!!   ***\n");
                        Msg_AsFem_Exit();
                        return false;
                    }
                    else{
                        elmtBlock._DomainName=str.substr(str.find_first_of("=")+1);
                        HasBlock=true;
                    }
                }
                else if(str.find("[end]")!=string::npos){
                    break;
                }
                else{
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: unknown option in [elmts] sub block                  !!!   ***\n");
                    Msg_AsFem_Exit();
                    return false;
                }
            }
            if(HasElmtBlock&&HasDofs&&HasElmtType){
                if(!HasMate){
                    elmtBlock._MateBlockName="";
                }
                if(!HasBlock) elmtBlock._DomainName="alldomain";
                elmtSystem.AddElmtBlock(elmtBlock);
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: [elmts] sub block information is not complete        !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        information is missing in [%35s] !!!   ***\n",elmtBlock._ElmtBlockName.c_str());
                if(!HasElmtType){
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type= found in [elmts] sub block                  !!!   ***\n");
                    return false;
                }

                if(!HasDofs){
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: no dofs= found in [elmts] sub block                  !!!   ***\n");
                    return false;
                }
                return false;
            }
        }
    }

    return true;
}