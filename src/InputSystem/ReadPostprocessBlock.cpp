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
//+++ Date   : 2021.02.22
//+++ Purpose: This function can read the [postprocess] block and its
//+++          subblock from our input file.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadPostprocessBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,Postprocess &postprocess,DofHandler &dofHandler){
    // postprocess block format
    // [postprocess]
    //  [pps1]
    //    type=ppstypename
    //    dof=u1 u2
    //    elmtid=mate1 [can be ignored]
    //    nodeid=all  [can be ignored]
    //    side=side1 side2
    //    domain=domain1 domain2
    //  [end]
    // [end]
    bool HasPPSBlock=false;
    bool HasPPSType=false;
    bool HasDof=false;
    bool HasElmtID=false;
    bool HasNodeID=false;

    PostprocessBlock ppsBlock;
    string tempstr,str0,substr;
    vector<double> number;
    vector<string> strlist;
    string msg;

    ppsBlock.Init();
    // now the str is [elmts]
    while(linenum<=lastendlinenum){
        getline(in,str);linenum+=1;
        str0=str;
        str=StringUtils::StrToLower(str0);
        str=StringUtils::RemoveStrSpace(str);
        if(StringUtils::IsCommentLine(str)||str.size()<1) continue;

        HasPPSBlock=false;
        HasPPSType=false;
        HasDof=false;


        if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
           (str.find("[end]")==string::npos&&str.find("[END]")==string::npos)){
            HasPPSBlock=false;
            tempstr=StringUtils::RemoveStrSpace(str0);
            if(tempstr.size()<3){
                MessagePrinter::PrintStars();
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no pps block name found in the [postprocess], a sub elmt block looks like [xxxx] should be given inside the [postprocess]",false);
                MessagePrinter::PrintStars();
                MessagePrinter::AsFem_Exit();
                return false;
            }
            else{
                ppsBlock._PPSBlockName=tempstr.substr(tempstr.find_first_of('[')+1,tempstr.find_first_of(']')-1);
                HasPPSBlock=true;
            }
            while(str.find("[end]")==string::npos&&str.find("[END]")==string::npos){
                getline(in,str);linenum+=1;
                str=StringUtils::RemoveStrSpace(str);
                str0=str;
                str=StringUtils::StrToLower(str0);

                if(StringUtils::IsCommentLine(str)||str.size()<1) {
                    //cout<<"str="<<str<<endl;
                    continue;
                }

                if(str.compare(0,5,"type=")==0){
                    substr=str.substr(str.find_first_of('=')+1);
                    if(substr.find("nodevalue")!=string::npos && substr.length()==9){
                        ppsBlock._PostprocessType=PostprocessType::NODALVALUEPPS;
                        ppsBlock._PostprocessTypeName="nodevalue";
                        HasPPSType=true;
                    }
                    else if(substr.find("elementvalue")!=string::npos && substr.length()==12){
                        ppsBlock._PostprocessType=PostprocessType::ELEMENTVALUEPPS;
                        ppsBlock._PostprocessTypeName="elementvalue";
                        HasPPSType=true;
                    }
                    else if(substr.find("elementalintegral")!=string::npos && substr.length()==17){
                        ppsBlock._PostprocessType=PostprocessType::ELEMENTINTEGRALPPS;
                        ppsBlock._PostprocessTypeName="elementalintegral";
                        HasPPSType=true;
                    }
                    else if(substr.find("sideintegral")!=string::npos && substr.length()==12){
                        ppsBlock._PostprocessType=PostprocessType::SIDEINTEGRALPPS;
                        ppsBlock._PostprocessTypeName="sideintegral";
                        HasPPSType=true;
                    }
                    else if(substr.find("area")!=string::npos && substr.length()==4){
                        ppsBlock._PostprocessType=PostprocessType::AREAPPS;
                        ppsBlock._PostprocessTypeName="area";
                        HasPPSType=true;
                    }
                    else if(substr.find("volume")!=string::npos && substr.length()==6){
                        ppsBlock._PostprocessType=PostprocessType::VOLUMEPPS;
                        ppsBlock._PostprocessTypeName="volume";
                        HasPPSType=true;
                    }
                    else if(substr.find("projvariablesideintegral")!=string::npos && substr.length()==24){
                        ppsBlock._PostprocessType=PostprocessType::PROJVARIABLESIDEINTEGRALPPS;
                        ppsBlock._PostprocessTypeName="projvariablesideintegral";
                        HasPPSType=true;
                    }
                    else if(substr.find("user")!=string::npos){
                        number=StringUtils::SplitStrNum(substr);
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        MessagePrinter::PrintErrorTxt("user postprocess is not supported yet in [postprocess] block, please waiting for the update",false);
                        MessagePrinter::PrintStars();
                        MessagePrinter::AsFem_Exit();
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
//                            int iuel=int(number[0]);
//                            elmtBlock._ElmtTypeName=str.substr(str.find_first_of('=')+1);
//                            if(iuel==1){
//                                elmtBlock._ElmtType=ElmtType::USER1ELMT;
//                            }
//                            else if(iuel==2){
//                                elmtBlock._ElmtType=ElmtType::USER2ELMT;
//                            }
//                            else if(iuel==3){
//                                elmtBlock._ElmtType=ElmtType::USER3ELMT;
//                            }
//                            else if(iuel==4){
//                                elmtBlock._ElmtType=ElmtType::USER4ELMT;
//                            }
//                            else if(iuel==5){
//                                elmtBlock._ElmtType=ElmtType::USER5ELMT;
//                            }
//                            else if(iuel==6){
//                                elmtBlock._ElmtType=ElmtType::USER6ELMT;
//                            }
//                            else if(iuel==7){
//                                elmtBlock._ElmtType=ElmtType::USER7ELMT;
//                            }
//                            else if(iuel==8){
//                                elmtBlock._ElmtType=ElmtType::USER8ELMT;
//                            }
//                            else if(iuel==9){
//                                elmtBlock._ElmtType=ElmtType::USER9ELMT;
//                            }
//                            else if(iuel==10){
//                                elmtBlock._ElmtType=ElmtType::USER10ELMT;
//                            }
//                            HasElmtType=true;
                        }
                    }
                    else{
                        msg="unsupported postprocess type in [postprocess] block, type="+substr+" is invalid";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                }
                else if(str.compare(0,4,"dof=")==0){
                    if(str.size()<4){
                        MessagePrinter::PrintStars();
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        MessagePrinter::PrintErrorTxt("no variable name found in [postprocess] sub block, 'variable=dof1 dof2 ...' is expected",false);
                        MessagePrinter::PrintStars();
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        substr=str0.substr(str0.find_first_of('=')+1);
                        string strtmp=substr;
                        if(!dofHandler.IsValidDofName(strtmp)){
                            MessagePrinter::PrintErrorTxt("dof name is not valid in [postprocess] sub block,the name should be the one list in the [dofs] block");
                            return false;
                        }
                        ppsBlock._VariableName=strtmp;
                        HasDof=true;
                    }
                }
                else if(str.compare(0,13,"projvariable=")==0){
                    if(str.size()<13){
                        MessagePrinter::PrintStars();
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        MessagePrinter::PrintErrorTxt("no projected variable name found in [postprocess] sub block, 'projvariable=variable-name' is expected",false);
                        MessagePrinter::PrintStars();
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        substr=str0.substr(str0.find_first_of('=')+1);
                        ppsBlock._ProjVariableName=substr;
                    }
                }
                else if(str.compare(0,11,"scalarmate=")==0){
                    if(str.size()<11){
                        MessagePrinter::PrintStars();
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        MessagePrinter::PrintErrorTxt("no scalar mate name found in [postprocess] sub block, 'scalarmate=mate-name' is expected",false);
                        MessagePrinter::PrintStars();
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        substr=str0.substr(str0.find_first_of('=')+1);
                        ppsBlock._ScalarMateName=substr;
                    }
                }
                else if(str.compare(0,11,"vectormate=")==0){
                    if(str.size()<11){
                        MessagePrinter::PrintStars();
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        MessagePrinter::PrintErrorTxt("no vector mate name found in [postprocess] sub block, 'vectormate=mate-name' is expected",false);
                        MessagePrinter::PrintStars();
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        substr=str0.substr(str0.find_first_of('=')+1);
                        ppsBlock._VectorMateName=substr;
                    }
                }
                else if(str.compare(0,10,"rank2mate=")==0){
                    if(str.size()<10){
                        MessagePrinter::PrintStars();
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        MessagePrinter::PrintErrorTxt("no rank-2 mate name found in [postprocess] sub block, 'rank2mate=mate-name' is expected",false);
                        MessagePrinter::PrintStars();
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        substr=str0.substr(str0.find_first_of('=')+1);
                        ppsBlock._VectorMateName=substr;
                    }
                }
                else if(str.compare(0,5,"side=")==0){
                    if(str.size()<5){
                        MessagePrinter::PrintStars();
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        MessagePrinter::PrintErrorTxt("no side name found in [postprocess] block, 'side=name_of_side' is expected",false);
                        MessagePrinter::PrintStars();
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        int i=str0.find_first_of('=');
                        substr=str0.substr(i+1,str0.length());
                        ppsBlock._BoundaryNameList=StringUtils::SplitStr(substr,' ');
                    }
                }
                else if(str.compare(0,7,"domain=")==0){
                    if(str.size()<7){
                        MessagePrinter::PrintStars();
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        MessagePrinter::PrintErrorTxt("no domain name found in [postprocess] block, 'domain=name_of_domain' is expected",false);
                        MessagePrinter::PrintStars();
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        int i=str0.find_first_of('=');
                        substr=str0.substr(i+1,str0.length());
                        ppsBlock._DomainNameList=StringUtils::SplitStr(substr,' ');
                    }
                }
                else if(str.compare(0,7,"elmtid=")==0){
                    if(!HasPPSType){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no 'type=' found in ["+ppsBlock._PPSBlockName+"] sub block, 'elmtid=' must be given after 'type=' in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    number=StringUtils::SplitStrNum(str);
                    if(number.size()<1){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no elmtid value found in ["+ppsBlock._PPSBlockName+"] sub block, 'elmtid=integer' should be given in [postprocess] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        ppsBlock._ElementID=static_cast<int>(number[0]);
                        HasElmtID=true;
                    }
                }
                else if(str.compare(0,7,"nodeid=")==0){
                    if(!HasPPSType){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no 'type=' found in ["+ppsBlock._PPSBlockName+"] sub block, 'nodeid=' must be given after 'type=' in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    number=StringUtils::SplitStrNum(str);
                    if(number.size()<1){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no nodeid value found in ["+ppsBlock._PPSBlockName+"] sub block, 'nodeid=integer' should be given in [postprocess] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        ppsBlock._NodeID=static_cast<int>(number[0]);
                        HasNodeID=true;
                    }
                }
                else if(str.compare(0,7,"iindex=")==0){
                    if(!HasPPSType){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no 'type=' found in ["+ppsBlock._PPSBlockName+"] sub block, 'iindex=' must be given after 'type=' in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    number=StringUtils::SplitStrNum(str);
                    if(number.size()<1){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no i-index value found in ["+ppsBlock._PPSBlockName+"] sub block, 'iindex=integer' should be given in [postprocess] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        ppsBlock._iInd=static_cast<int>(number[0]);
                        if(ppsBlock._iInd>3||ppsBlock._iInd<1){
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            msg="i-index="+to_string(ppsBlock._iInd)+" is invalid";
                            MessagePrinter::PrintErrorTxt(msg);
                            MessagePrinter::AsFem_Exit();
                        }
                    }
                }
                else if(str.compare(0,7,"jindex=")==0){
                    if(!HasPPSType){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no 'type=' found in ["+ppsBlock._PPSBlockName+"] sub block, 'jindex=' must be given after 'type=' in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    number=StringUtils::SplitStrNum(str);
                    if(number.size()<1){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no j-index value found in ["+ppsBlock._PPSBlockName+"] sub block, 'jindex=integer' should be given in [postprocess] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        ppsBlock._jInd=static_cast<int>(number[0]);
                        if(ppsBlock._jInd>3||ppsBlock._jInd<1){
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            msg="j-index="+to_string(ppsBlock._jInd)+" is invalid";
                            MessagePrinter::PrintErrorTxt(msg);
                            MessagePrinter::AsFem_Exit();
                        }
                    }
                }
                else if(str.compare(0,10,"component=")==0){
                    if(!HasPPSType){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no 'type=' found in ["+ppsBlock._PPSBlockName+"] sub block, 'component=' must be given after 'type=' in [bcs] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                    }
                    number=StringUtils::SplitStrNum(str);
                    if(number.size()<1){
                        MessagePrinter::PrintErrorInLineNumber(linenum);
                        msg="no component value found in ["+ppsBlock._PPSBlockName+"] sub block, 'component=integer' should be given in [postprocess] sub block";
                        MessagePrinter::PrintErrorTxt(msg);
                        MessagePrinter::AsFem_Exit();
                        return false;
                    }
                    else{
                        ppsBlock._Component=static_cast<int>(number[0]);
                        if(ppsBlock._Component>3||ppsBlock._Component<1){
                            MessagePrinter::PrintErrorInLineNumber(linenum);
                            msg="component="+to_string(ppsBlock._Component)+" is invalid";
                            MessagePrinter::PrintErrorTxt(msg);
                            MessagePrinter::AsFem_Exit();
                        }
                    }
                }
                else if(str.find("[end]")!=string::npos){
                    break;
                }
                else{
                    MessagePrinter::PrintStars();
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("unknown options in [postprocess] sub block",false);
                    MessagePrinter::PrintStars();
                    MessagePrinter::AsFem_Exit();
                    return false;
                }
            }
            if(HasPPSBlock&&HasPPSType){
                postprocess.AddPostprocessBlock(ppsBlock);
                ppsBlock.Init();
                if((!HasDof)||(!HasElmtID)||(!HasNodeID)){}
            }
            else{
                msg="information in [postprocess] sub block is not complete, information is missing in ["+ppsBlock._PPSBlockName+"]";
                MessagePrinter::PrintErrorTxt(msg,false);
                if(!HasPPSType){
                    msg="no 'type=' option can be found in ["+ppsBlock._PPSBlockName+"] sub block";
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