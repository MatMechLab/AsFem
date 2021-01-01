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
//+++ Date   : 2020.07.12
//+++ Purpose: This function can read the [nonlinearsolver] block 
//+++          from our input file.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadNonlinearSolverBlock(ifstream &in,string str,int &linenum,NonlinearSolver &nonlinearSolver){
    // dof block format:
    // [linearsolver]
    //   type=lu [gmres]
    //   maxiters=10000
    //   tol=1.0e-9
    // [end]
    

    bool HasType=false;
    vector<double> numbers;
    string namelist;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;

    

    _nonlinearSolverBlock.Init();// use the default value
    
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(StringUtils::IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            str=StringUtils::StrToLower(str);
            continue;
        }
        if(str.find("type=")!=string::npos||
           str.find("TYPE=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=StringUtils::RemoveStrSpace(substr);
            if((substr.find("nr")!=string::npos||substr.find("NR")!=string::npos)&&
               substr.length()==2){
                HasType=true;
                _nonlinearSolverBlock._SolverTypeName="newton-raphson";
                _nonlinearSolverBlock._SolverType=NonlinearSolverType::NEWTON;
            }
            else if((substr.find("newtonls")!=string::npos||substr.find("NEWTONLS")!=string::npos)&&
               substr.length()==8){
                HasType=true;
                _nonlinearSolverBlock._SolverTypeName="newton with line search";
                _nonlinearSolverBlock._SolverType=NonlinearSolverType::NEWTONLS;
            }
            else if((substr.find("newtontr")!=string::npos||substr.find("NEWTONTR")!=string::npos)&&
               substr.length()==8){
                HasType=true;
                _nonlinearSolverBlock._SolverTypeName="newton with trust region";
                _nonlinearSolverBlock._SolverType=NonlinearSolverType::NEWTONTR;
            }
            else if((substr.find("bfgs")!=string::npos||substr.find("BFGS")!=string::npos)&&
               substr.length()==4){
                HasType=true;
                _nonlinearSolverBlock._SolverTypeName="BFGS";
                _nonlinearSolverBlock._SolverType=NonlinearSolverType::BFGS;
            }
            else if((substr.find("broyden")!=string::npos||substr.find("BROYDEN")!=string::npos)&&
               substr.length()==7){
                HasType=true;
                _nonlinearSolverBlock._SolverTypeName="broyden";
                _nonlinearSolverBlock._SolverType=NonlinearSolverType::BROYDEN;
            }
            else if((substr.find("ngmres")!=string::npos||substr.find("NGMRES")!=string::npos)&&
               substr.length()==6){
                HasType=true;
                _nonlinearSolverBlock._SolverTypeName="ngmres";
                _nonlinearSolverBlock._SolverType=NonlinearSolverType::NEWTONGMRES;
            }
            else if((substr.find("ncg")!=string::npos||substr.find("NCG")!=string::npos)&&
               substr.length()==6){
                HasType=true;
                _nonlinearSolverBlock._SolverTypeName="ncg";
                _nonlinearSolverBlock._SolverType=NonlinearSolverType::NEWTONCG;
            }
            else{
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("unsupported solver type in the [nonlinearsolver] block");
                MessagePrinter::AsFem_Exit();
            }
        }
        else if(str.find("maxiters")!=string::npos||str.find("MAXITERS")!=string::npos){
            if(!HasType){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no 'type=' is found in the [nonlinearsolver] block, 'maxiters=' must be given after 'type='");
                MessagePrinter::AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=StringUtils::SplitStrNum(substr);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("maxiters= number can not be found in the [nonlinearsolver] block, maxiters=integer is expected");
                MessagePrinter::AsFem_Exit();
            }
            else{
                if(int(numbers[0])<1){
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("invalid maxiters number in the [nonlinearsolver] block, maxiters=integer is expected",false);
                    MessagePrinter::AsFem_Exit();
                }
                _nonlinearSolverBlock._MaxIters=int(numbers[0]);
            }
        }
        else if(str.find("r_rel_tol=")!=string::npos||
                 str.find("R_rel_tol=")!=string::npos||
                 str.find("R_Rel_Tol=")!=string::npos||
                 str.find("R_REL_TOL=")!=string::npos){
            if(!HasType){
                MessagePrinter::PrintErrorTxt("no 'type=' found in the [nonlinearsolver] block, r_rel_tol= should be given after 'type='");
                MessagePrinter::AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=StringUtils::SplitStrNum(substr);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no r_rel_tol found in the [nonlinearsolver] block, r_rel_tol=real should be given");
                MessagePrinter::AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("invalid r_rel_tol found in [nonlinearsolver] block, r_rel_tol=real should be given");
                    MessagePrinter::AsFem_Exit();
                }
                _nonlinearSolverBlock._RRelTol=numbers[0];
            }
        }
        else if(str.find("r_abs_tol=")!=string::npos||
                 str.find("R_abs_tol=")!=string::npos||
                 str.find("R_Abs_Tol=")!=string::npos||
                 str.find("R_ABS_TOL=")!=string::npos){
            if(!HasType){
                MessagePrinter::PrintErrorTxt("no 'type=' found in the [nonlinearsolver] block, 'r_abs_tol=' should be given after 'type='",false);
                MessagePrinter::AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=StringUtils::SplitStrNum(substr);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no 'r_abs_tol' found in [nonlinearsolver] block, 'r_abs_tol=' should be given after 'type='",false);
                MessagePrinter::AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("invalid r_abs_tol found in [nonlinearsolver] block, 'r_abs_tol=real' should be given",false);
                    MessagePrinter::AsFem_Exit();
                }
                _nonlinearSolverBlock._RAbsTol=numbers[0];
            }
        }
        else if(str.find("stol=")!=string::npos||
                 str.find("Stol=")!=string::npos||
                 str.find("STOL=")!=string::npos){
            if(!HasType){
                MessagePrinter::PrintErrorTxt("no 'type=' found in [nonlinearsolver] block, 'stol=' should be given after 'type='",false);
                MessagePrinter::AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=StringUtils::SplitStrNum(substr);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no stol= number found in [nonlinearsolver] block, stol=real should be given",false);
                MessagePrinter::AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("invalid stol= number found in [nonlinearsolver] block, stol=real should be given",false);
                    MessagePrinter::AsFem_Exit();
                }
                _nonlinearSolverBlock._STol=numbers[0];
            }
        }
        else if(str.find("[]")!=string::npos){
            MessagePrinter::PrintErrorInLineNumber(linenum);
            MessagePrinter::PrintErrorTxt("the bracket pair is not complete in the [nonlinearsolver] block",false);
            MessagePrinter::AsFem_Exit();
        }
        else{
            MessagePrinter::PrintErrorInLineNumber(linenum);
            MessagePrinter::PrintErrorTxt("unknown option in [nonlinearsolver] block",false);
            MessagePrinter::AsFem_Exit();
        }
        getline(in,str);linenum+=1;
    }

    nonlinearSolver.SetOptionsFromNonlinearSolverBlock(_nonlinearSolverBlock);
    
    return HasType;
}