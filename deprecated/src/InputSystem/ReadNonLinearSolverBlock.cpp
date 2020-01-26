#include "InputSystem/InputSystem.h"

bool InputSystem::ReadNonLinearSolverBlock(ifstream &in,string str,int &linenum,NonLinearSolverBlockInfo &nonlinearSolverBlockInfo){
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

    nonlinearSolverBlockInfo.Reset();// use the default value
    
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            continue;
        }
        if(str.find("type=")!=string::npos||
           str.find("TYPE=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if((substr.find("nr")!=string::npos||substr.find("NR")!=string::npos)&&
               substr.length()==2){
                HasType=true;
                nonlinearSolverBlockInfo._SolverType=NonLinearSolverType::NEWTONRAPHSON;
            }
            else if((substr.find("newton")!=string::npos||substr.find("NEWTON")!=string::npos)&&
               substr.length()==6){
                HasType=true;
                nonlinearSolverBlockInfo._SolverType=NonLinearSolverType::NEWTONRAPHSON;
            }
            else if((substr.find("newtonraphson")!=string::npos||
                     substr.find("NewtonRaphson")!=string::npos||
                     substr.find("NEWTONRAPHSON")!=string::npos)&&
                   substr.length()==13){
                HasType=true;
                nonlinearSolverBlockInfo._SolverType=NonLinearSolverType::NEWTONRAPHSON;
            }
            else if((substr.find("linear")!=string::npos||
                     substr.find("Linear")!=string::npos||
                     substr.find("Linear")!=string::npos)&&
                   substr.length()==6){
                HasType=true;
                nonlinearSolverBlockInfo._SolverType=NonLinearSolverType::LINEAR;
            }
            else if((substr.find("nrlinesearch")!=string::npos||
                     substr.find("NrLineSearch")!=string::npos||
                     substr.find("NRLINESEARCH")!=string::npos)&&
                   substr.length()==12){
                HasType=true;
                nonlinearSolverBlockInfo._SolverType=NonLinearSolverType::NEWTONWITHLINESEARCH;
            }
            else if((substr.find("nrbisectlinesearch")!=string::npos||
                     substr.find("NrBisectLineSearch")!=string::npos||
                     substr.find("NRBISECTLINESEARCH")!=string::npos)&&
                   substr.length()==18){
                HasType=true;
                nonlinearSolverBlockInfo._SolverType=NonLinearSolverType::NEWTONWITHBISECTIONLINESEARCH;
            }
            else if((substr.find("nrregulalinesearch")!=string::npos||
                     substr.find("NrRegulaLineSearch")!=string::npos||
                     substr.find("NRREGULALINESEARCH")!=string::npos)&&
                   substr.length()==18){
                HasType=true;
                nonlinearSolverBlockInfo._SolverType=NonLinearSolverType::NEWTONWITHREGULARFALSILINESEARCH;
            }
            else if((substr.find("nrmodify")!=string::npos||
                     substr.find("NrModify")!=string::npos||
                     substr.find("NRMODIFY")!=string::npos)&&
                    substr.length()==8){
                HasType=true;
                nonlinearSolverBlockInfo._SolverType=NonLinearSolverType::MODIFIEDNEWTON;
            }
            else{
                Msg_Input_LineError(linenum);
                Msg_Input_NonLinearSolverBlock_InvalidSolverOption();
                Msg_AsFem_Exit();
            }
        }
        else if(str.find("maxiters")!=string::npos||str.find("MAXITERS")!=string::npos){
            if(!HasType){
                cout<<"*** Error: no type= is found in [nonlinearsolver] block !!!***"<<endl;
                cout<<"***        maxiters should be given after type  !!!        ***"<<endl;
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_NonLinearSolverBlock_MaxIterNotFound();
                Msg_AsFem_Exit();
            }
            else{
                if(int(numbers[0])<1){
                    Msg_Input_LineError(linenum);
                    Msg_Input_NonLinearSolverBlock_MaxIterInvalid();
                    Msg_AsFem_Exit();
                }
                nonlinearSolverBlockInfo._MaxIters=int(numbers[0]);
            }
        }
        else if(str.find("maxsearch")!=string::npos||
                str.find("MaxSearch")!=string::npos||
                str.find("MAXSEARCH")!=string::npos){
            if(!HasType){
                cout<<"*** Error: no type= is found in [nonlinearsolver] block !!!***"<<endl;
                cout<<"***        maxsearch should be given after type  !!!       ***"<<endl;
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                cout<<"*** Error: maxsearch not found in [nonlinearsolver] block!!***"<<endl;
                cout<<"***        maxsearch=integer should be given after type!!! ***"<<endl;
                Msg_AsFem_Exit();
            }
            else{
                if(int(numbers[0])<1){
                    Msg_Input_LineError(linenum);
                    Msg_Input_NonLinearSolverBlock_MaxIterInvalid();
                    cout<<"*** Error: maxsearch= invalid in [nonlinearsolver] block!!!***"<<endl;
                    cout<<"***        maxsearch=integer should be given after type!!! ***"<<endl;
                    Msg_AsFem_Exit();
                }
                nonlinearSolverBlockInfo._MaxLineSearch=int(numbers[0]);
            }
        }
        else if(str.find("r_rel_tol=")!=string::npos||
                 str.find("R_rel_tol=")!=string::npos||
                 str.find("R_Rel_Tol=")!=string::npos||
                 str.find("R_REL_TOL=")!=string::npos){
            if(!HasType){
                cout<<"*** Error: no type= is found in [linearsolver] block  !!!  ***"<<endl;
                cout<<"***        r_rel_tol should be given after type !!!        ***"<<endl;
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_NonLinearSolverBlock_RrelTolNotFound();
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    Msg_Input_LineError(linenum);
                    Msg_Input_NonLinearSolverBlock_RrelTolInvalid();
                    Msg_AsFem_Exit();
                }
                nonlinearSolverBlockInfo._Rreltol=numbers[0];
            }
        }
        else if(str.find("r_abs_tol")!=string::npos||
                 str.find("R_Abs_Tol")!=string::npos||
                 str.find("R_abs_tol")!=string::npos||
                 str.find("R_ABS_TOL")!=string::npos){
            if(!HasType){
                cout<<"*** Error: no type= is found in [linearsolver] block  !!!  ***"<<endl;
                cout<<"***        r_abs_tol should be given after type !!!        ***"<<endl;
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_NonLinearSolverBlock_RabsTolNotFound();
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    Msg_Input_LineError(linenum);
                    Msg_Input_NonLinearSolverBlock_RabsTolInvalid();
                    Msg_AsFem_Exit();
                }
                nonlinearSolverBlockInfo._Rabstol=numbers[0];
            }
        }
        else if(str.find("e_rel_tol")!=string::npos||
                str.find("E_rel_tol")!=string::npos||
                str.find("E_Rel_Tol")!=string::npos||
                str.find("E_REL_TOL")!=string::npos){
            if(!HasType){
                cout<<"*** Error: no type= is found in [linearsolver] block  !!!  ***"<<endl;
                cout<<"***        e_rel_tol should be given after type !!!        ***"<<endl;
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_NonLinearSolverBlock_ErelTolNotFound();
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    Msg_Input_LineError(linenum);
                    Msg_Input_NonLinearSolverBlock_ErelTolInvalid();
                    Msg_AsFem_Exit();
                }
                nonlinearSolverBlockInfo._Ereltol=numbers[0];
            }
        }
        else if(str.find("e_abs_tol")!=string::npos||
                str.find("E_Abs_Tol")!=string::npos||
                str.find("E_ABS_TOL")!=string::npos){
            if(!HasType){
                cout<<"*** Error: no type= is found in [linearsolver] block  !!!  ***"<<endl;
                cout<<"***        e_abs_tol should be given after type !!!        ***"<<endl;
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_NonLinearSolverBlock_EabsTolNotFound();
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    Msg_Input_LineError(linenum);
                    Msg_Input_NonLinearSolverBlock_EabsTolInvalid();
                    Msg_AsFem_Exit();
                }
                nonlinearSolverBlockInfo._Eabstol=numbers[0];
            }
        }
        else if(str.find("stol=")!=string::npos||
                str.find("STol=")!=string::npos||
                str.find("STOL=")!=string::npos){
            if(!HasType){
                cout<<"*** Error: no type= is found in [linearsolver] block  !!!  ***"<<endl;
                cout<<"***        stol should be given after type !!!             ***"<<endl;
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_NonLinearSolverBlock_STolNotFound();
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    Msg_Input_LineError(linenum);
                    Msg_Input_NonLinearSolverBlock_STolInvalid();
                    Msg_AsFem_Exit();
                }
                nonlinearSolverBlockInfo._Stol=numbers[0];
            }
        }
        else if(str.find("[]")!=string::npos){
            Msg_Input_LineError(linenum);
            Msg_Input_BlockBracketNotComplete();
            Msg_AsFem_Exit();
        }
        else{
            Msg_Input_LineError(linenum);
            cout<<"*** Error: unknown option in [nonlinearsolver] block!!!    ***"<<endl;
            Msg_AsFem_Exit();
        }
        getline(in,str);linenum+=1;

    }
    return HasType;
}