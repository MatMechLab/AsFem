//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadMateBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,MateSystem &mateSystem){
    // each bc block should looks like:
    //   [mate1]
    //     type=elastic
    //     params=1.0 2.0
    //   [end]
    // important: now , str already contains [mate] !!!

    bool HasMateBlock=false;
    MateBlock mateBlock;
    string tempstr,str0;
    vector<double> number;

    bool HasElmt=false;
    bool HasValue=false;
    bool HasBlock=false;

    mateBlock.Reset();

    // now str="[mate]"
    while (linenum<=lastendlinenum){
        getline(in,str);linenum+=1;
        str=StrToLower(str);
        str=RemoveStrSpace(str);
        if(IsCommentLine(str)||str.size()<1) continue;

        

        if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
          (str.find("[end]")==string::npos&&str.find("[END]")==string::npos)){
            tempstr=RemoveStrSpace(str);
            mateBlock.Reset();
            HasElmt=false;
            HasValue=false;
            HasBlock=false;
            if(tempstr.size()<3){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no material block name found in [mates] sub block    !!!   ***\n");
                return false;
            }
            else{
                mateBlock._MateBlockName=tempstr.substr(1,tempstr.size()-1-1);
                HasBlock=true;
            }
            while(str.find("[end]")==string::npos&&str.find("[END]")==string::npos){
                getline(in,str);linenum+=1;
                str=StrToLower(str);
                str0=str;
                str=RemoveStrSpace(str0);
                if(IsCommentLine(str)||str.size()<1) continue;
                if(str.find("type=")!=string::npos){
                    string substr=str.substr(str.find_first_of('=')+1);
                    if(substr.find("constpoisson")!=string::npos&&substr.length()==12){
                        mateBlock._MateTypeName="constpoisson";
                        mateBlock._MateType=MateType::ConstPoissonMate;
                        HasElmt=true;
                    }
                    else if(substr.find("nlpoisson")!=string::npos&&substr.length()==9){
                        mateBlock._MateTypeName="nlpoisson";
                        mateBlock._MateType=MateType::NonLinearPoissonMate;
                        HasElmt=true;
                    }
                    else if(substr.find("tensorpoisson")!=string::npos&&substr.length()==13){
                        mateBlock._MateTypeName="tensorpoisson";
                        mateBlock._MateType=MateType::TensorPoissonMate;
                        HasElmt=true;
                    }
                    else if(substr.find("constdiffusion")!=string::npos&&substr.length()==14){
                        mateBlock._MateTypeName="constdiffusion";
                        mateBlock._MateType=MateType::ConstDiffusionMate;
                        HasElmt=true;
                    }
                    else if(substr.find("nldiffusion")!=string::npos&&substr.length()==11){
                        mateBlock._MateTypeName="nldiffusion";
                        mateBlock._MateType=MateType::NonLinearDiffusionMate;
                        HasElmt=true;
                    }
                    else if(substr.find("cahnhilliard")!=string::npos&&substr.length()==12){
                        mateBlock._MateTypeName="cahnhilliard";
                        mateBlock._MateType=MateType::CahnHilliardMate;
                        HasElmt=true;
                    }
                    else if(substr.find("tensorcahnhilliard")!=string::npos&&substr.length()==18){
                        mateBlock._MateTypeName="tensorcahnhilliard";
                        mateBlock._MateType=MateType::TensorCahnHilliardMate;
                        HasElmt=true;
                    }
                    else if(substr.find("linearelastic")!=string::npos&&substr.length()==13){
                        mateBlock._MateTypeName="linearelastic";
                        mateBlock._MateType=MateType::LinearElasticMate;
                        HasElmt=true;
                    }
                    else if(substr.find("saintvenant")!=string::npos&&substr.length()==11){
                        mateBlock._MateTypeName="saintvenant";
                        mateBlock._MateType=MateType::SaintVenantMate;
                        HasElmt=true;
                    }
                    else if(substr.find("neohookean")!=string::npos&&substr.length()==10){
                        mateBlock._MateTypeName="neohookean";
                        mateBlock._MateType=MateType::NeoHookeanMate;
                        HasElmt=true;
                    }
                    else if(substr.find("mooneyrivlin")!=string::npos&&substr.length()==12){
                        mateBlock._MateTypeName="mooneyrivlin";
                        mateBlock._MateType=MateType::MooneyRivlinMate;
                        HasElmt=true;
                    }
                    else if(substr.find("constwave")!=string::npos&&substr.length()==9){
                        mateBlock._MateTypeName="constwave";
                        mateBlock._MateType=MateType::ConstWaveMate;
                        HasElmt=true;
                    }
                    else if(substr.find("thermelastic")!=string::npos&&substr.length()==12){
                        mateBlock._MateTypeName="thermelastic";
                        mateBlock._MateType=MateType::LinearThermalMechanicsMate;
                        HasElmt=true;
                    }
                    else if(substr.find("lpffracture")!=string::npos&&substr.length()==11){
                        mateBlock._MateTypeName="lpffracture";
                        mateBlock._MateType=MateType::LinearElasticPhaseFieldFracMate;
                        HasElmt=true;
                    }
                    else if(substr.find("cohesivepffrac")!=string::npos&&substr.length()==14){
                        mateBlock._MateTypeName="cohesivepffrac";
                        mateBlock._MateType=MateType::CohesivePFFracMate;
                        HasElmt=true;
                    }
                    else if(substr.find("anisopffrac")!=string::npos&&substr.length()==11){
                        mateBlock._MateTypeName="anisopffrac";
                        mateBlock._MateType=MateType::AnisoLinearElasticPhaseFieldFracMate;
                        HasElmt=true;
                    }
                    else if(substr.find("miehelinear")!=string::npos&&substr.length()==11){
                        mateBlock._MateTypeName="miehelinear";
                        mateBlock._MateType=MateType::MieheLinearElasticFracMate;
                        HasElmt=true;
                    }
                    else if(substr.find("mieheneohookean")!=string::npos&&substr.length()==15){
                        mateBlock._MateTypeName="mieheneohookean";
                        mateBlock._MateType=MateType::MieheNeoHookeanFracMate;
                        HasElmt=true;
                    }
                    else if(substr.find("bordenlinear")!=string::npos&&substr.length()==12){
                        mateBlock._MateTypeName="bordenlinear";
                        mateBlock._MateType=MateType::BordenLinearElasticFracMate;
                        HasElmt=true;
                    }
                    else if(substr.find("elasticch")!=string::npos&&substr.length()==9){
                        mateBlock._MateTypeName="elasticch";
                        mateBlock._MateType=MateType::LinearElasticCahnHilliardMate;
                        HasElmt=true;
                    }
                    else if(substr.find("dendrite")!=string::npos&&substr.length()==8){
                        mateBlock._MateTypeName="dendrite";
                        mateBlock._MateType=MateType::DendriteMate;
                        HasElmt=true;
                    }
                    else if(substr.find("currentthermal")!=string::npos&&substr.length()==14){
                        mateBlock._MateTypeName="currentthermal";
                        mateBlock._MateType=MateType::CurrentThermalMate;
                        HasElmt=true;
                    }
                    else if(substr.find("user")!=string::npos){
                        number=SplitStrNum(substr);
                        if(number.size()<1){
                            Msg_Input_LineError(linenum);
                            PetscPrintf(PETSC_COMM_WORLD,"*** Error: no user number found umat in [mates] sub block       !!!   ***\n");
                            PetscPrintf(PETSC_COMM_WORLD,"***        type=user1,user2,...,user20 is expected              !!!   ***\n");
                            Msg_AsFem_Exit();
                            return false;
                        }
                        
                        if(int(number[0])<1||int(number[0])>20){
                            Msg_Input_LineError(linenum);
                            PetscPrintf(PETSC_COMM_WORLD,"*** Error: user number is invalid in [mates] sub block          !!!   ***\n");
                            PetscPrintf(PETSC_COMM_WORLD,"***        type=user1,user2,...,user20 is expected              !!!   ***\n");
                            Msg_AsFem_Exit();
                            return false;
                        }
                        else{
                            mateBlock._MateTypeName=str.substr(str.find_first_of('=')+1);
                            switch(int(number[0])){
                                case 1:
                                  mateBlock._MateType=MateType::User1Mate;
                                  break;
                                case 2:
                                  mateBlock._MateType=MateType::User2Mate;
                                  break;
                                case 3:
                                  mateBlock._MateType=MateType::User3Mate;
                                  break;
                                case 4:
                                  mateBlock._MateType=MateType::User4Mate;
                                  break;
                                case 5:
                                  mateBlock._MateType=MateType::User5Mate;
                                  break;
                                case 6:
                                  mateBlock._MateType=MateType::User6Mate;
                                  break;
                                case 7:
                                  mateBlock._MateType=MateType::User7Mate;
                                  break;
                                case 8:
                                  mateBlock._MateType=MateType::User8Mate;
                                  break;
                                case 9:
                                  mateBlock._MateType=MateType::User9Mate;
                                  break;
                                case 10:
                                  mateBlock._MateType=MateType::User10Mate;
                                  break;
                                case 11:
                                  mateBlock._MateType=MateType::User11Mate;
                                  break;
                                case 12:
                                  mateBlock._MateType=MateType::User12Mate;
                                  break;
                                case 13:
                                  mateBlock._MateType=MateType::User13Mate;
                                  break;
                                case 14:
                                  mateBlock._MateType=MateType::User14Mate;
                                  break;
                                case 15:
                                  mateBlock._MateType=MateType::User15Mate;
                                  break;
                                case 16:
                                  mateBlock._MateType=MateType::User16Mate;
                                  break;
                                case 17:
                                  mateBlock._MateType=MateType::User17Mate;
                                  break;
                                case 18:
                                  mateBlock._MateType=MateType::User18Mate;
                                  break;
                                case 19:
                                  mateBlock._MateType=MateType::User19Mate;
                                  break;
                                case 20:
                                  mateBlock._MateType=MateType::User20Mate;
                                  break;
                                default:
                                  mateBlock._MateType=MateType::NullMate;
                                  break;
                            }
                            HasElmt=true;
                        }
                    }
                    else{
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported material in [mates] sub block            !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        type= %-25s is invalid       !!!   ***\n",substr.c_str());
                        Msg_AsFem_Exit();
                    }
                }
                else if(str0.find("params=")!=string::npos){
                    if(!HasElmt){
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type= found in [mates] sub block                  !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        params must be given after type in [mates] sub block !!!   ***\n");
                        Msg_AsFem_Exit();
                    }
                    number=SplitStrNum(str0);
                    if(number.size()<1){
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: parameters not found in [mates] sub block            !!!   ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        params=x y z... is expected                          !!!   ***\n");
                        return false;
                    }
                    else{
                        mateBlock._Params.clear();
                        mateBlock._Params=number;
                        HasValue=true;
                    }
                }
            }
            if(HasBlock&&HasElmt){
                HasBlock=true;
                if(!HasValue) {
                    mateBlock._Params.clear();
                    mateBlock._Params.push_back(1.0);
                }
                mateSystem.AddMateBlock(mateBlock);
                HasMateBlock=true;
            }
            else{
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: information is not complete in [mates] sub block     !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        information is missing in [%-20s] !!!   ***\n",mateBlock._MateBlockName.c_str());
                if(!HasElmt){
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type found in [mates] sub block                   !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        type= must be given in each [mates] sub block        !!!   ***\n");
                }
                return false;
            }
        }
    }

    return HasMateBlock;
}