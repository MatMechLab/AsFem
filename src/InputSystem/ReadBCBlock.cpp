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
//+++ Date   : 2020.07.10
//+++ Purpose: This function can read the [bcs] block and its 
//+++          subblock from our input file.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
    string msg;

    bool HasElmt=false;
    bool HasValue=false;
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
            tempstr=StringUtils::RemoveStrSpace(str);
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
                if(str.find("type=")!=string::npos){
                    string substr=str.substr(str.find_first_of('=')+1);
                    if(substr.find("dirichlet")!=string::npos && substr.length()==9){
                        bcblock._BCTypeName="dirichlet";
                        bcblock._BCType=BCType::DIRICHLETBC;
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
                    else if(substr.find("user")!=string::npos){
                        number=StringUtils::SplitStrNum(str);
                        if(number.size()<1){
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            msg="no user number found umat in ["+bcblock._BCBlockName+"] sub block, 'type=user1,user2,...,user10' is expected";
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
                else if(str.find("dof=")!=string::npos){
                    if(!HasElmt){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="type= is not given yet in ["+bcblock._BCBlockName+"] sub block, 'dof=' should be given after 'type=' in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    if(str.size()<5){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="dof name not found in ["+bcblock._BCBlockName+"] sub block, 'dof=' should be given after 'type=' in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    else{
                        bcblock._DofName=str.substr(4);
                        if(!dofHandler.IsValidDofName(bcblock._DofName)){
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            msg="invalid dof name found in ["+bcblock._BCBlockName+"] sub block, dof must be one of the names list in [dofs] block";
                            MessagePrinter::PrintErrorTxt(msg);
                            MessagePrinter::AsFem_Exit();
                            return false;
                        }
                        bcblock._DofID=dofHandler.GetDofIDviaDofName(bcblock._DofName);
                        
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
                        if(bclist.size()>0)
                        {
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
            }
            if(HasDof&&HasBoundary&&HasElmt&&HasBCBlock){
                HasBCBlock=true;
                if(!HasValue) bcblock._BCValue=0.0;
                bcSystem.AddBCBlock2List(bcblock);
            }
            else{
                msg="information is not complete in [mates] sub block, some information is missing in ["+bcblock._BCBlockName+"]";
                MessagePrinter::PrintErrorTxt(msg,false);
                if(!HasElmt){
                    msg="no type found in ["+bcblock._BCBlockName+"] sub block, 'type=' must be given in each [bcs] sub block";
                    MessagePrinter::PrintErrorTxt(msg);
                }

                if(!HasDof){
                    msg="no dof found in ["+bcblock._BCBlockName+"] sub block, 'dof=dof_name' must be given in each [bcs] sub block";
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