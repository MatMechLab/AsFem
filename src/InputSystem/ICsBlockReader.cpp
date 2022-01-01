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
//+++ Date   : 2021.08.05
//+++ Purpose: Implement the reader for [ics] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/ICsBlockReader.h"


void ICsBlockReader::PrintHelper(){
    MessagePrinter::PrintStars(MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("The complete information for [ics] block:",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("[ics]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  [ic-1]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    type=ic-type-1",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    dof=dof-name-1",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    params=val1 val2",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    domain=domain-name-1",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  [end]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  [ic-2]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    type=ic-type-2",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    dof=dof-name-2",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    params=val3 val4",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("    domain=domain-name-2",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  [end]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("[end]",MessageColor::BLUE);
    MessagePrinter::PrintStars(MessageColor::BLUE);
}
//**********************************************************
bool ICsBlockReader::ReadICBlock(ifstream &in, string str, const int &lastendlinenum, int &linenum, ICSystem &icSystem, DofHandler &dofHandler){
    bool HasICBlock=false;
    ICBlock icblock;
    string tempstr,str0,stro;
    vector<double> number;
    vector<string> domainlist;
    string msg;

    bool HasElmt=false;
    bool HasValue=false;
    bool HasDof=false;
    bool HasDomain=false;

    // now str="[bcs]"
    while (linenum<=lastendlinenum){
        getline(in,str);linenum+=1;
        str0=str;
        str=StringUtils::RemoveStrSpace(str);
        stro=str;
        str=StringUtils::StrToLower(str);
        if(StringUtils::IsCommentLine(str)||str.size()<1) continue;

        HasElmt=false;
        HasValue=false;
        HasDof=false;

        icblock.Init();

        if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
          (str.find("[end]")==string::npos&&str.find("[END]")==string::npos)){
            HasICBlock=false;
            str=StringUtils::StrToLower(str);
            tempstr=StringUtils::RemoveStrSpace(str);
            if(tempstr.size()<3){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                msg="no IC block name found in the [ics] sub block";
                MessagePrinter::PrintErrorTxt(msg,false);
                MessagePrinter::AsFem_Exit();
                return false;
            }
            else{
                icblock.Init();
                icblock._ICBlockName=tempstr.substr(1,tempstr.size()-1-1);
                HasICBlock=true;
            }
            while(str.find("[end]")==string::npos&&str.find("[END]")==string::npos){
                getline(in,str);linenum+=1;
                str0=str;
                str=StringUtils::RemoveStrSpace(str);
                stro=str;
                str=StringUtils::StrToLower(stro);
                if(StringUtils::IsCommentLine(str)||str.size()<1) continue;
                
                if(str.find("helper")!=string::npos){
                    PrintHelper();
                    return false;
                }
                else if(str.find("type=")!=string::npos){
                    string substr=str.substr(str.find('=')+1);
                    if(substr.find("const")!=string::npos&&substr.length()==5){
                        icblock._ICTypeName="const";
                        icblock._ICType=ICType::CONSTIC;
                        HasElmt=true;
                    }
                    else if(substr.find("random")!=string::npos&&substr.length()==6){
                        icblock._ICTypeName="random";
                        icblock._ICType=ICType::RANDOMIC;
                        HasElmt=true;
                    }
                    else if(substr.find("circle")!=string::npos&&substr.length()==6){
                        icblock._ICTypeName="circle";
                        icblock._ICType=ICType::CIRCLEIC;
                        HasElmt=true;
                    }
                    else if(substr.find("smoothcircle")!=string::npos&&substr.length()==12){
                        icblock._ICTypeName="smoothcircle";
                        icblock._ICType=ICType::SMOOTHCIRCLEIC;
                        HasElmt=true;
                    }
                    else if(substr.find("sphere")!=string::npos&&substr.length()==6){
                        icblock._ICTypeName="sphere";
                        icblock._ICType=ICType::SPHERICALIC;
                        HasElmt=true;
                    }
                    else if(substr.find("rectangle")!=string::npos&&substr.length()==9){
                        icblock._ICTypeName="rectangle";
                        icblock._ICType=ICType::RECTANGLEIC;
                        HasElmt=true;
                    }
                    else if(substr.find("cubic")!=string::npos&&substr.length()==5){
                        icblock._ICTypeName="cubic";
                        icblock._ICType=ICType::CUBICIC;
                        HasElmt=true;
                    }
                    else if(substr.find("user")!=string::npos){
                        number=StringUtils::SplitStrNum(str);
                        if(number.size()<1){
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            msg="no user number found in 'type=user' in the ["+icblock._ICBlockName+"] sub block, 'type=user1,user2,...,user10' is expected";
                            MessagePrinter::PrintErrorTxt(msg,false);
                            MessagePrinter::AsFem_Exit();
                            return false;
                        }

                        if(int(number[0])<1||int(number[0])>10){
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            msg="user number is invalid in the ["+icblock._ICBlockName+"] sub block, 'type=user1,user2,...,user10' is expected";
                            MessagePrinter::PrintErrorTxt(msg,false);
                            MessagePrinter::AsFem_Exit();
                            return false;
                        }
                        else{
                            icblock._ICTypeName=str.substr(str.find_first_of('=')+1);
                            switch(int(number[0])){
                                case 1:
                                    icblock._ICType=ICType::USER1IC;
                                    break;
                                case 2:
                                    icblock._ICType=ICType::USER2IC;
                                    break;
                                case 3:
                                    icblock._ICType=ICType::USER3IC;
                                    break;
                                case 4:
                                    icblock._ICType=ICType::USER4IC;
                                    break;
                                case 5:
                                    icblock._ICType=ICType::USER5IC;
                                    break;
                                case 6:
                                    icblock._ICType=ICType::USER6IC;
                                    break;
                                case 7:
                                    icblock._ICType=ICType::USER7IC;
                                    break;
                                case 8:
                                    icblock._ICType=ICType::USER8IC;
                                    break;
                                case 9:
                                    icblock._ICType=ICType::USER9IC;
                                    break;
                                case 10:
                                    icblock._ICType=ICType::USER10IC;
                                    break;
                                default:
                                    icblock._ICType=ICType::NULLIC;
                                    break;
                            }
                            HasElmt=true;
                        }
                    }
                    else{
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="unsupported IC type in the ["+icblock._ICBlockName+"] sub block, type="+substr+" is invalid";
                        MessagePrinter::PrintErrorTxt(msg,false);
                        MessagePrinter::AsFem_Exit();
                    }
                }
                else if(str.find("dof=")!=string::npos){
                    if(!HasElmt){
                        msg="'type=' is not given yet in ["+icblock._ICBlockName+"] sub block, 'dof=' should be given after 'type=' in [ics] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    if(str.size()<5){
                        msg="no dof name found in the ["+icblock._ICBlockName+"] sub block, 'dof=dof_name' should be given in [ics] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    else{
                        icblock._DofName=stro.substr(4);
                        if(!dofHandler.IsValidDofName(icblock._DofName)){
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            msg="invalid dof name found in ["+icblock._ICBlockName+"] sub block, it must be one of the names list in the [dofs] block";
                            MessagePrinter::PrintErrorTxt(msg,false);
                            MessagePrinter::AsFem_Exit();

                            return false;
                        }
                        icblock._DofID=dofHandler.GetDofIDviaDofName(icblock._DofName);
                        HasDof=true;
                    }
                }
                else if(str.find("domain=")!=string::npos){
                    if(!HasElmt){
                        msg="'type=' is not given yet in ["+icblock._ICBlockName+"] sub block, 'domoain=' should be given after 'type=' in [ics] sub block";
                        MessagePrinter::PrintErrorTxt(msg,false);
                        MessagePrinter::AsFem_Exit();
                    }
                    if(str.size()<7){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no domain name found in ["+icblock._ICBlockName+"] sub block, 'domain=domain_name' should be given in [ics] sub block";
                        MessagePrinter::PrintErrorTxt(msg,false);
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        int i=str0.find_first_of('=');
                        string substr=str0.substr(i+1,str0.length());
                        domainlist=StringUtils::SplitStr(substr,' ');
                        if(domainlist.size()>0){
                            if(StringUtils::IsUniqueStrVec(domainlist)){
                                icblock._DomainNameList=domainlist;
                                HasDomain=true;
                            }
                            else{
                                MessagePrinter::PrintErrorInLineNumber(linenum);
                                msg="duplicated domain name in the ["+icblock._ICBlockName+"] sub block";
                                MessagePrinter::PrintErrorTxt(msg,false);
                                MessagePrinter::AsFem_Exit();
                            }
                        }
                        else{
                            HasDomain=false;
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            msg="no domain name found in ["+icblock._ICBlockName+"] sub block, 'domain=domain1 domain2 ...' is expected";
                            MessagePrinter::PrintErrorTxt(msg);
                            MessagePrinter::AsFem_Exit();
                        }
                    }
                }
                else if(str.find("params=")!=string::npos){
                    if(!HasElmt){
                        msg="type= is not given yet in ["+icblock._ICBlockName+"] sub block, 'params=' should be given after 'type=' in [ics] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    number=StringUtils::SplitStrNum(str0);
                    if(number.size()<1){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no ic parameters are found in ["+icblock._ICBlockName+"] sub block, 'params=real' value should be given in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg,false);
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        icblock._Parameters=number;
                        HasValue=true;
                    }
                }
            }
            if(HasDof&&HasElmt&&HasICBlock){
                HasICBlock=true;
                if(!HasValue) {
                    icblock._Parameters.clear();
                    icblock._Parameters.push_back(0.0);
                    icblock._Parameters.push_back(1.0);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Warning: no params are given in [ics] sub block             !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***          default value will be used in your IC element      !!!   ***\n");
                    msg="no params are given in ["+icblock._ICBlockName+"] sub block, default value will be used in your IC element";
                    MessagePrinter::PrintWarningTxt(msg);
                }
                if(!HasDomain){
                    icblock._DomainNameList.clear();
                    icblock._DomainNameList.push_back("alldomain");
                    HasDomain=true;
                }
                icSystem.AddICBlock2List(icblock);
                icblock.Init();
            }
            else{
                msg="information is not complete in ["+icblock._ICBlockName+"] sub block, some information is missing, please check your input file";
                MessagePrinter::PrintErrorTxt(msg,false);
                if(!HasElmt){
                    msg="no type found in ["+icblock._ICBlockName+"] sub block, 'type=' must be given in each [ics] sub block";
                    MessagePrinter::PrintErrorTxt(msg,false);
                    MessagePrinter::AsFem_Exit();
                }

                if(!HasDof){
                    msg="no dof name found in ["+icblock._ICBlockName+"] sub block, ' dof=dof-name' must be given in each [ics] sub block";
                    MessagePrinter::PrintErrorTxt(msg,false);
                    MessagePrinter::AsFem_Exit();
                }
                return false;
            }
        }
    }

    return true;
}
