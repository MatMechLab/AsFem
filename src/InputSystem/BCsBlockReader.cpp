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
//+++ Date   : 2021.08.05
//+++ Purpose: Implement the reader for [bcs] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/BCsBlockReader.h"

#include "DofHandler/DofHandler.h"

void BCsBlockReader::PrintHelper(){
    MessagePrinter::PrintStars(MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("The complete information for [bcs] block:",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("[bcs]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  [bc-1]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    type=bc-type-1",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    dofs=dof1 dofs2 ...",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    value=bc-value-1",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    params=val1 val2 ...",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    boundary=boundary-name-1",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  [end]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  [bc-2]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    type=bc-type-2",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    dofs=dof2",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    value=bc-value-2",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    params=v1 v2 ...",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    boundary=boundary-name-2",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  [end]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("[end]",MessageColor::BLUE);
    MessagePrinter::PrintStars(MessageColor::BLUE);
}

//******************************************************
bool BCsBlockReader::ReadBCBlock(ifstream &in, string str, const int &lastendlinenum, int &linenum, BCSystem &bcSystem, DofHandler &dofHandler){
    bool HasBCBlock=false;
    BCBlock bcblock;
    string tempstr,str0,substr;
    vector<double> number;
    vector<string> bclist;
    string msg;

    bool HasElmt=false;
    bool HasValue=false;
    bool HasParams=false;
    bool HasBoundary=false;
    bool HasDof=false;

    bcblock.Init();

    // now str="[bcs]"
    while (linenum<=lastendlinenum){
        getline(in,str);linenum+=1;
        str0=str;
        str=StringUtils::StrToLower(str0);
        str=StringUtils::RemoveStrSpace(str);
        if(StringUtils::IsCommentLine(str)||str.size()<1) continue;

        HasElmt=false;
        HasValue=false;
        HasBoundary=false;
        HasDof=false;

        if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
          (str.find("[end]")==string::npos&&str.find("[END]")==string::npos)){
            HasBCBlock=false;
            tempstr=StringUtils::RemoveStrSpace(str0);
            if(tempstr.size()<3){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no boundary block name can be found in [bcs] sub block");
                MessagePrinter::AsFem_Exit();
                return false;
            }
            else{
                bcblock.Init();
                bcblock._BCBlockName=tempstr.substr(1,tempstr.size()-1-1);
                HasBCBlock=true;
            }
            while(str.find("[end]")==string::npos&&str.find("[END]")==string::npos){
                getline(in,str);linenum+=1;
                str0=str;
                str=StringUtils::RemoveStrSpace(str0);
                str=StringUtils::StrToLower(str);
                if(StringUtils::IsCommentLine(str)||str.size()<1) continue;

                if(str.find("helper")!=string::npos){
                    PrintHelper();
                    return false;
                }
                else if(str.find("type=")!=string::npos){
                    string substr=str.substr(str.find_first_of('=')+1);
                    if(substr.find("dirichlet")!=string::npos && substr.length()==9){
                        bcblock._BCTypeName="dirichlet";
                        bcblock._BCType=BCType::DIRICHLETBC;
                        HasElmt=true;
                    }
                    else if(substr.find("user1dirichlet")!=string::npos && substr.length()==14){
                        bcblock._BCTypeName="user1dirichlet";
                        bcblock._BCType=BCType::USER1DIRICHLETBC;
                        HasElmt=true;
                    }
                    else if(substr.find("user2dirichlet")!=string::npos && substr.length()==14){
                        bcblock._BCTypeName="user2dirichlet";
                        bcblock._BCType=BCType::USER2DIRICHLETBC;
                        HasElmt=true;
                    }
                    else if(substr.find("user3dirichlet")!=string::npos && substr.length()==14){
                        bcblock._BCTypeName="user3dirichlet";
                        bcblock._BCType=BCType::USER3DIRICHLETBC;
                        HasElmt=true;
                    }
                    else if(substr.find("user4dirichlet")!=string::npos && substr.length()==14){
                        bcblock._BCTypeName="user4dirichlet";
                        bcblock._BCType=BCType::USER4DIRICHLETBC;
                        HasElmt=true;
                    }
                    else if(substr.find("user5dirichlet")!=string::npos && substr.length()==14){
                        bcblock._BCTypeName="user5dirichlet";
                        bcblock._BCType=BCType::USER5DIRICHLETBC;
                        HasElmt=true;
                    }
                    else if(substr.find("neumann")!=string::npos && substr.length()==7){
                        bcblock._BCTypeName="neumann";
                        bcblock._BCType=BCType::NEUMANNBC;
                        HasElmt=true;
                    }
                    else if(substr.find("nodaldirichlet")!=string::npos && substr.length()==14){
                        bcblock._BCTypeName="nodaldirichlet";
                        bcblock._BCType=BCType::NODALDIRICHLETBC;
                        HasElmt=true;
                    }
                    else if(substr.find("nodalneumann")!=string::npos && substr.length()==12){
                        bcblock._BCTypeName="nodalneumann";
                        bcblock._BCType=BCType::NODALNEUMANNBC;
                        HasElmt=true;
                    }
                    else if(substr.find("pressure")!=string::npos && substr.length()==8){
                        bcblock._BCTypeName="pressure";
                        bcblock._BCType=BCType::PRESSUREBC;
                        HasElmt=true;
                    }
                    else if(substr.find("user")!=string::npos && substr.size()<=6){
                        number=StringUtils::SplitStrNum(str);
                        if(number.size()<1){
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            msg="no user number found in ["+bcblock._BCBlockName+"] sub block, 'type=user1,user2,...,user10' is expected";
                            MessagePrinter::PrintErrorTxt(msg);
                            MessagePrinter::AsFem_Exit();
                            return false;
                        }

                        if(int(number[0])<1||int(number[0])>10){
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            msg="user number is invalid in ["+bcblock._BCBlockName+"] sub block, 'type=user1,user2,...,user10' is expected";
                            MessagePrinter::PrintErrorTxt(msg);
                            MessagePrinter::AsFem_Exit();

                            return false;
                        }
                        else{
                            bcblock._BCTypeName=str.substr(str.find_first_of('=')+1);
                            switch(int(number[0])){
                                case 1:
                                  bcblock._BCType=BCType::USER1BC;
                                  break;
                                case 2:
                                  bcblock._BCType=BCType::USER2BC;
                                  break;
                                case 3:
                                  bcblock._BCType=BCType::USER3BC;
                                  break;
                                case 4:
                                  bcblock._BCType=BCType::USER4BC;
                                  break;
                                case 5:
                                  bcblock._BCType=BCType::USER5BC;
                                  break;
                                case 6:
                                  bcblock._BCType=BCType::USER6BC;
                                  break;
                                case 7:
                                  bcblock._BCType=BCType::USER7BC;
                                  break;
                                case 8:
                                  bcblock._BCType=BCType::USER8BC;
                                  break;
                                case 9:
                                  bcblock._BCType=BCType::USER9BC;
                                  break;
                                case 10:
                                  bcblock._BCType=BCType::USER10BC;
                                  break;
                                default:
                                  bcblock._BCType=BCType::NULLBC;
                                  break;
                            }
                            HasElmt=true;
                        }
                    }
                    else{
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="unsupported type in ["+bcblock._BCBlockName+"] sub block, 'type="+substr+"' is invalid";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                }
                else if(str.find("dofs=")!=string::npos){
                    if(!HasElmt){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="type= is not given yet in ["+bcblock._BCBlockName+"] sub block, 'dofs=' should be given after 'type=' in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    if(str.size()<6){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="dof name not found in ["+bcblock._BCBlockName+"] sub block, 'dofs=' should be given after 'type=' in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    else{
                        int i=str0.find_first_of('=');
                        substr=str0.substr(i+1,str0.length());
                        bcblock._DofsName=StringUtils::SplitStr(substr,' ');
                        bcblock._DofIDs.resize(bcblock._DofsName.size());
                        i=0;
                        for(const auto &name:bcblock._DofsName){
                            if(!dofHandler.IsValidDofName(name)){
                                MessagePrinter::PrintErrorInLineNumber(linenum);
                                msg="invalid dof name found in ["+bcblock._BCBlockName+"] sub block, dof must be one of the names list in [dofs] block";
                                MessagePrinter::PrintErrorTxt(msg);
                                MessagePrinter::AsFem_Exit();
                                return false;
                            }
                            bcblock._DofIDs[i]=dofHandler.GetDofIDviaDofName(name);
                            i+=1;
                        }
                        HasDof=true;
                    }
                }
                else if(str.find("boundary=")!=string::npos){
                    if(!HasElmt){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="'type=' is not given yet in ["+bcblock._BCBlockName+"] sub block, 'boundary=' should be given after 'type=' in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    if(str.size()<10){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="boundary name not found in ["+bcblock._BCBlockName+"] sub block, 'boundary=bcname' should be given in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        int i=str0.find_first_of('=');
                        string substr=str0.substr(i+1,str0.length());
                        bclist=StringUtils::SplitStr(substr,' ');
                        if(bclist.size()>0){
                            if(StringUtils::IsUniqueStrVec(bclist)){
                                bcblock._BoundaryNameList=bclist;
                                HasBoundary=true;
                            }
                            else{
                                MessagePrinter::PrintErrorInLineNumber(linenum);
                                msg="duplicated boundary name in ["+bcblock._BCBlockName+"] sub block";
                                MessagePrinter::PrintErrorTxt(msg);
                                MessagePrinter::AsFem_Exit();
                            }
                        }
                        else{
                            HasBoundary=false;
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            msg="no boundary name found in ["+bcblock._BCBlockName+"] sub block, 'boundary=bc1 bc2 ...' is expected";
                            MessagePrinter::PrintErrorTxt(msg);
                            MessagePrinter::AsFem_Exit();
                        }
                    }
                }
                else if(str.find("value=")!=string::npos){
                    if(!HasElmt){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no 'type=' found in ["+bcblock._BCBlockName+"] sub block, 'value=' must be given after 'type=' in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    number=StringUtils::SplitStrNum(str);
                    if(number.size()<1&&str.find("t")==string::npos){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no bc value found in ["+bcblock._BCBlockName+"] sub block, 'value=real value' should be given in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        HasValue=true;
                        int i;
                        i=str.find_first_of("=");
                        string substr=str.substr(i+1);
                        if(substr.find("t")!=string::npos&&substr.find("*")==string::npos){
                            bcblock._IsTimeDependent=true;
                            bcblock._BCValue=1.0;
                        }
                        else if(substr.find("-t")!=string::npos&&substr.find("*")==string::npos){
                            bcblock._IsTimeDependent=true;
                            bcblock._BCValue=-1.0;
                        }
                        else{
                            if(StringUtils::IsValidExpression(substr)){
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
                else if(str.find("params=")!=string::npos){
                    if(!HasElmt){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no 'type=' found in ["+bcblock._BCBlockName+"] sub block, 'value=' must be given after 'type=' in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    number=StringUtils::SplitStrNum(str0);
                    if(number.size()<1){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="parameters can not be found in ["+bcblock._BCBlockName+"] block, 'params=x y z..' is expected";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        bcblock._Parameters=number;
                        HasParams=true;
                    }
                }
                else if(str.find("[end]")!=string::npos){
                    break;
                }
                else{
                    MessagePrinter::PrintStars();
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("unknown options in ["+bcblock._BCBlockName+"] sub block",false);
                    MessagePrinter::PrintStars();
                    MessagePrinter::AsFem_Exit();
                    return false;
                }
            }
            if(HasDof&&HasBoundary&&HasElmt&&HasBCBlock){
                HasBCBlock=true;
                if(!HasValue) bcblock._BCValue=0.0;
                if(!HasParams) bcblock._Parameters.push_back(1.0);
                bcSystem.AddBCBlock2List(bcblock);
                bcblock.Init();
            }
            else{
                msg="information is not complete in [bcs] sub block, some information is missing in ["+bcblock._BCBlockName+"]";
                MessagePrinter::PrintErrorTxt(msg,false);
                if(!HasElmt){
                    msg="no type found in ["+bcblock._BCBlockName+"] sub block, 'type=' must be given in each [bcs] sub block";
                    MessagePrinter::PrintErrorTxt(msg);
                }

                if(!HasDof){
                    msg="no dof found in ["+bcblock._BCBlockName+"] sub block, 'dofs=dof_name' must be given in each [bcs] sub block";
                    MessagePrinter::PrintErrorTxt(msg);
                }

                if(!HasBoundary){
                    msg="no boundary found in ["+bcblock._BCBlockName+"] sub block, 'boundary=bcname' must be given in each [bcs] sub block";
                    MessagePrinter::PrintErrorTxt(msg);
                }
                return false;
            }
        }
    }

    return true;
}
