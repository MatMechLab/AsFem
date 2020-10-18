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
//+++ Date   : 2020.07.01
//+++ Purpose: This function can read the [elmts] block and its 
//+++          subblock from our input file.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
    string msg;


    elmtBlock.Init();
    // now the str is [elmts]
    while(linenum<=lastendlinenum){
        getline(in,str);linenum+=1;
        str0=str;
        str=StringUtils::StrToLower(str0);
        str=StringUtils::RemoveStrSpace(str);
        if(StringUtils::IsCommentLine(str)||str.size()<1) continue;

        HasElmtBlock=false;
        HasElmtType=false;
        HasDofs=false;
        HasBlock=false;
        HasMate=false;

        if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
          (str.find("[end]")==string::npos&&str.find("[END]")==string::npos)){
            HasElmtBlock=false;
            tempstr=StringUtils::RemoveStrSpace(str);
            if(tempstr.size()<3){
                MessagePrinter::PrintStars();
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no elmt block name found in the [elmts], a sub elmt block looks like [xxxx] should be given inside the [elmts]",false);
                MessagePrinter::PrintStars();
                MessagePrinter::AsFem_Exit();
                return false;
            }
            else{
                elmtBlock._ElmtBlockName=tempstr.substr(tempstr.find_first_of('[')+1,tempstr.find_first_of(']')-1);
                HasElmtBlock=true;
            }
            while(str.find("[end]")==string::npos&&str.find("[END]")==string::npos){
                getline(in,str);linenum+=1;
                str0=str;
                str=StringUtils::StrToLower(str0);
                str=StringUtils::RemoveStrSpace(str);

                if(StringUtils::IsCommentLine(str)||str.size()<1) {
                    //cout<<"str="<<str<<endl;
                    continue;
                }

                if(str.find("type=")!=string::npos){
                    substr=str.substr(str.find_first_of('=')+1);
                    if(substr.find("poisson")!=string::npos && substr.length()==7){
                        elmtBlock._ElmtTypeName="poisson";
                        elmtBlock._ElmtType=ElmtType::POISSONELMT;
                        HasElmtType=true;
                    }
                    else if(substr.find("diffusion")!=string::npos && substr.length()==9){
                        elmtBlock._ElmtTypeName="diffusion";
                        elmtBlock._ElmtType=ElmtType::DIFFUSIONELMT;
                        HasElmtType=true;
                    }
                    else if(substr.find("wave")!=string::npos && substr.length()==4){
                        elmtBlock._ElmtTypeName="wave";
                        elmtBlock._ElmtType=ElmtType::WAVEELMT;
                        HasElmtType=true;
                    }
                    else if(substr.find("mechanics")!=string::npos && substr.length()==9){
                        elmtBlock._ElmtTypeName="mechanics";
                        elmtBlock._ElmtType=ElmtType::MECHANICSELMT;
                        HasElmtType=true;
                    }
                    else if(substr.find("cahnhilliard")!=string::npos && substr.length()==12){
                        elmtBlock._ElmtTypeName="cahnhilliard";
                        elmtBlock._ElmtType=ElmtType::CAHNHILLIARDELMT;
                        HasElmtType=true;
                    }
                    else if(substr.find("mechcahnhilliard")!=string::npos && substr.length()==16){
                        elmtBlock._ElmtTypeName="mechcahnhilliard";
                        elmtBlock._ElmtType=ElmtType::MECHCAHNHILLIARDELMT;
                        HasElmtType=true;
                    }
                    else if(substr.find("thermalmechanics")!=string::npos && substr.length()==16){
                        elmtBlock._ElmtTypeName="thermalmechanics";
                        elmtBlock._ElmtType=ElmtType::THERMALMECHANICSELMT;
                        HasElmtType=true;
                    }
                    else if(substr.find("miehefrac")!=string::npos && substr.length()==9){
                        elmtBlock._ElmtTypeName="miehefrac";
                        elmtBlock._ElmtType=ElmtType::MIEHEFRACELMT;
                        HasElmtType=true;
                    }
                    else if(substr.find("thermal")!=string::npos && substr.length()==7){
                        elmtBlock._ElmtTypeName="thermal";
                        elmtBlock._ElmtType=ElmtType::THERMALCONDUCTELMT;
                        HasElmtType=true;
                    }
                    else if(substr.find("dendrite")!=string::npos && substr.length()==8){
                        elmtBlock._ElmtTypeName="dendrite";
                        elmtBlock._ElmtType=ElmtType::DENDRITEELMT;
                        HasElmtType=true;
                    }
                    else if(substr.find("user")!=string::npos){
                        number=StringUtils::SplitStrNum(substr);
                        if(number.size()<1){
                            MessagePrinter::PrintStars();
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            MessagePrinter::PrintErrorTxt("no user element number can be found in the [elmts] block, 'type=user1,user2,user3'... is expected",false);
                            MessagePrinter::PrintStars();
                            MessagePrinter::AsFem_Exit();

                            return false;
                        }
                        if(int(number[0])<1||int(number[0])>20){
                            MessagePrinter::PrintStars();
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            MessagePrinter::PrintErrorTxt("invalid user element number in [elmts] block, type=user1~user20  is expected",false);
                            MessagePrinter::PrintStars();
                            MessagePrinter::AsFem_Exit();
                            return false;
                        }
                        else{
                            int iuel=int(number[0]);
                            elmtBlock._ElmtTypeName=str.substr(str.find_first_of('=')+1);
                            if(iuel==1){
                                elmtBlock._ElmtType=ElmtType::USER1ELMT;
                            }
                            else if(iuel==2){
                                elmtBlock._ElmtType=ElmtType::USER2ELMT;
                            }
                            else if(iuel==3){
                                elmtBlock._ElmtType=ElmtType::USER3ELMT;
                            }
                            else if(iuel==4){
                                elmtBlock._ElmtType=ElmtType::USER4ELMT;
                            }
                            else if(iuel==5){
                                elmtBlock._ElmtType=ElmtType::USER5ELMT;
                            }
                            else if(iuel==6){
                                elmtBlock._ElmtType=ElmtType::USER6ELMT;
                            }
                            else if(iuel==7){
                                elmtBlock._ElmtType=ElmtType::USER7ELMT;
                            }
                            else if(iuel==8){
                                elmtBlock._ElmtType=ElmtType::USER8ELMT;
                            }
                            else if(iuel==9){
                                elmtBlock._ElmtType=ElmtType::USER9ELMT;
                            }
                            else if(iuel==10){
                                elmtBlock._ElmtType=ElmtType::USER10ELMT;
                            }
                            else if(iuel==11){
                                elmtBlock._ElmtType=ElmtType::USER11ELMT;
                            }
                            else if(iuel==12){
                                elmtBlock._ElmtType=ElmtType::USER12ELMT;
                            }
                            else if(iuel==13){
                                elmtBlock._ElmtType=ElmtType::USER13ELMT;
                            }
                            else if(iuel==14){
                                elmtBlock._ElmtType=ElmtType::USER14ELMT;
                            }
                            else if(iuel==15){
                                elmtBlock._ElmtType=ElmtType::USER15ELMT;
                            }
                            else if(iuel==16){
                                elmtBlock._ElmtType=ElmtType::USER16ELMT;
                            }
                            else if(iuel==17){
                                elmtBlock._ElmtType=ElmtType::USER17ELMT;
                            }
                            else if(iuel==18){
                                elmtBlock._ElmtType=ElmtType::USER18ELMT;
                            }
                            else if(iuel==19){
                                elmtBlock._ElmtType=ElmtType::USER19ELMT;
                            }
                            else if(iuel==20){
                                elmtBlock._ElmtType=ElmtType::USER20ELMT;
                            }
                            HasElmtType=true;
                        }
                    }
                    else{
                        msg="unsupported element type in [elmts] block, type="+substr+" is invalid";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                }
                else if(str.find("dofs=")!=string::npos){
                    if(str.size()<5){
                        MessagePrinter::PrintStars();
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        MessagePrinter::PrintErrorTxt("no dof name found in [elmts] sub block, 'dofs=dof1 dof2 ...' is expected",false);
                        MessagePrinter::PrintStars();
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        substr=str0.substr(str0.find_first_of('=')+1);
                        vector<string> strtmp=StringUtils::SplitStr(substr,' ');
                        if(!dofHandler.IsValidDofNameVec(strtmp)){
                            MessagePrinter::PrintErrorTxt("dof name is not valid in [elmts] sub block,the name should be the one list in the [dofs] block");
                            return false;
                        }
                        if(StringUtils::IsUniqueStrVec(strtmp)){
                            elmtBlock._DofsNameList=StringUtils::SplitStr(substr,' ');
                            elmtBlock._DofsIDList=dofHandler.GetDofsIndexFromNameVec(elmtBlock._DofsNameList);
                            elmtBlock._nDofs=static_cast<int>(elmtBlock._DofsIDList.size());
                            HasDofs=true;
                        }
                        else{
                            HasDofs=false;
                            MessagePrinter::PrintStars();
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            MessagePrinter::PrintErrorTxt("duplicated dof name found in [elmts] sub block, 'dofs=dof1 dof2 ...' is expected",false);
                            MessagePrinter::PrintStars();
                            MessagePrinter::AsFem_Exit();
                            return false;
                        }
                    }
                }
                else if(str.find("mate=")!=string::npos){
                    if(str.size()<5){
                        MessagePrinter::PrintStars();
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        MessagePrinter::PrintErrorTxt("material block name not found in [elmts] sub block, 'mate=name_of_material_block' is expected",false);
                        MessagePrinter::PrintStars();
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        elmtBlock._MateBlockName=str.substr(str.find_first_of("=")+1);
                        HasMate=true;
                    }
                }
                else if(str.find("domain=")!=string::npos){
                    if(str.size()<7){
                        MessagePrinter::PrintStars();
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        MessagePrinter::PrintErrorTxt("no domain name found in [elmts] block, 'domain=name_of_domain' is expected",false);
                        MessagePrinter::PrintStars();
                        MessagePrinter::AsFem_Exit();
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
                    MessagePrinter::PrintStars();
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("unknown options in [elmts] sub block",false);
                    MessagePrinter::PrintStars();
                    MessagePrinter::AsFem_Exit();
                    return false;
                }
            }
            if(HasElmtBlock&&HasDofs&&HasElmtType){
                if(!HasMate){
                    elmtBlock._MateBlockName="";
                }
                if(!HasBlock) elmtBlock._DomainName="alldomain";
                elmtSystem.AddBulkElmtBlock2List(elmtBlock);
            }
            else{
                msg="information in [elmts] sub block is not complete, information is missing in ["+elmtBlock._ElmtBlockName+"]";
                MessagePrinter::PrintErrorTxt(msg,false);
                if(!HasElmtType){
                    msg="no 'type=' option can be found in ["+elmtBlock._ElmtBlockName+"] sub block";
                    MessagePrinter::PrintErrorTxt(msg,false);
                    MessagePrinter::AsFem_Exit();
                    return false;
                }

                if(!HasDofs){
                    msg="no 'dofs=' option can be found in ["+elmtBlock._ElmtBlockName+"] sub block";
                    MessagePrinter::PrintErrorTxt(msg,false);
                    MessagePrinter::AsFem_Exit();
                    return false;
                }
                return false;
            }
        }
    }

    return true;
}