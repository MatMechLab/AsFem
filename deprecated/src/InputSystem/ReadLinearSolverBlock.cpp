#include "InputSystem/InputSystem.h"

bool InputSystem::ReadLinearSolverBlock(ifstream &in,string str,int &linenum,LinearSolverBlockInfo &linearSolverBlockInfo)
{
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

    linearSolverBlockInfo.Reset();// use the default value
    
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos)
    {
        if(IsCommentLine(str)||str.length()<1)
        {
            getline(in,str);linenum+=1;
            continue;
        }
        if(str.find("type=")!=string::npos||
           str.find("TYPE=")!=string::npos)
        {
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if((substr.find("lu")!=string::npos||substr.find("LU")!=string::npos)&&
               substr.length()==2){
                HasType=true;
                linearSolverBlockInfo._SolverName="lu";
                linearSolverBlockInfo._SolverType=LinearSolverType::SPARSELU;
            }
            else if((substr.find("qr")!=string::npos||
                     substr.find("QR")!=string::npos)&&
                    substr.length()==2){
                HasType=true;
                linearSolverBlockInfo._SolverName="qr";
                linearSolverBlockInfo._SolverType=LinearSolverType::SPARSEQR;
            }
            else if((substr.find("pardiso")!=string::npos||
                     substr.find("PARDISO")!=string::npos)&&
                    substr.length()==7){
                HasType=true;
                linearSolverBlockInfo._SolverName="pardiso";
                linearSolverBlockInfo._SolverType=LinearSolverType::PARDISOLU;
                #ifndef USE_PARDISO
                cout<<"*** Warning: you try to use pardiso, but it isn\'t enabled! ***"<<endl;
                cout<<"***          make sure you enabled pardiso in CMakeLists!  ***"<<endl;
                cout<<"***          besides, you must install intel icc compiler!!***"<<endl;
                cout<<"***          AsFem choose SparseLU for you, instead!!!     ***"<<endl;
                linearSolverBlockInfo._SolverType=LinearSolverType::SPARSELU;
                #endif
            }
            else if((substr.find("umfpack")!=string::npos||
                     substr.find("Umfpack")!=string::npos||
                     substr.find("UmfPack")!=string::npos||
                     substr.find("UMFPACK")!=string::npos)&&
                    substr.length()==7){
                HasType=true;
                linearSolverBlockInfo._SolverName="umfpack";
                linearSolverBlockInfo._SolverType=LinearSolverType::UMFPACK;
                #ifndef USE_SUITESPARSE
                cout<<"*** Warning: you try to use umfpack, but it isn\'t enabled! ***"<<endl;
                cout<<"***          make sure you enabled pardiso in CMakeLists!  ***"<<endl;
                cout<<"***          besides, you must install intel icc compiler!!***"<<endl;
                cout<<"***          AsFem choose SparseLU for you, instead!!!     ***"<<endl;
                linearSolverBlockInfo._SolverType=LinearSolverType::SPARSELU;
                #endif
            }
            else if((substr.find("cg")!=string::npos||
                     substr.find("CG")!=string::npos)&&
                    substr.length()==2){
                HasType=true;
                linearSolverBlockInfo._SolverName="cg";
                linearSolverBlockInfo._SolverType=LinearSolverType::CG;
            }
            else if((substr.find("bicg")!=string::npos||
                     substr.find("BICG")!=string::npos)&&
                    substr.length()==4){
                HasType=true;
                linearSolverBlockInfo._SolverName="bicg";
                linearSolverBlockInfo._SolverType=LinearSolverType::BICG;
            }
            else if((substr.find("gmres")!=string::npos||
                     substr.find("GMRES")!=string::npos)&&
                    substr.length()==5){
                HasType=true;
                linearSolverBlockInfo._SolverName="gmres";
                linearSolverBlockInfo._SolverType=LinearSolverType::GMRES;
            }
            else
            {
                Msg_Input_LineError(linenum);
                Msg_Input_LinearSolverBlock_InvalidSolverOption();
                HasType=false;
                Msg_AsFem_Exit();
            }
            
        }
        else if(str.find("maxiters")!=string::npos||str.find("MAXITERS")!=string::npos){
            if(!HasType){
                cout<<"*** Error: no type= is found in [linearsolver] block  !!!  ***"<<endl;
                cout<<"***        maxiters should be given after type=  !!!       ***"<<endl;
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_LinearSolverBlock_MaxIterNotFound();
                Msg_AsFem_Exit();
            }
            else{
                if(int(numbers[0])<1){
                    Msg_Input_LineError(linenum);
                    Msg_Input_LinearSolverBlock_MaxIterInvalid();
                    Msg_AsFem_Exit();
                }
                linearSolverBlockInfo._MaxIters=int(numbers[0]);
            }
        }
        else if(str.find("tol")!=string::npos||str.find("TOL")!=string::npos){
            if(!HasType){
                cout<<"*** Error: no type= is found in [linearsolver] block  !!!  ***"<<endl;
                cout<<"***        tol should be given after type !!!              ***"<<endl;
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_LinearSolverBlock_TolNotFound();
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<1.0e-14){
                    Msg_Input_LineError(linenum);
                    Msg_Input_LinearSolverBlock_TolInvalid();
                    Msg_AsFem_Exit();
                }
                linearSolverBlockInfo._Tol=numbers[0];
            }
        }
        else if(str.find("restart")!=string::npos||
                str.find("Restart")!=string::npos||
                str.find("RESTART")!=string::npos){
            if(!HasType){
                cout<<"*** Error: no type= is found in [linearsolver] block  !!!  ***"<<endl;
                cout<<"***        restart should be given after type !!!          ***"<<endl;
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                cout<<"*** Error: restart number not found in [linearsolver] block***"<<endl;
                cout<<"***        restart=positive integer is expected !!!        ***"<<endl;
                Msg_AsFem_Exit();
            }
            else{
                if(int(numbers[0])<10){
                    Msg_Input_LineError(linenum);
                    cout<<"*** Error: invalid restart number in [linearsolver] block!!***"<<endl;
                    cout<<"***        restart=positive integer is expected !!!        ***"<<endl;
                    Msg_AsFem_Exit();
                }
                linearSolverBlockInfo._Restart=int(numbers[0]);
            }
        }
        else if(str.find("[]")!=string::npos)
        {
            Msg_Input_LineError(linenum);
            Msg_Input_BlockBracketNotComplete();
            Msg_AsFem_Exit();
        }
        else{
            Msg_Input_LineError(linenum);
            cout<<"*** Error: unknown option in [linearsolver] block     !!!  ***"<<endl;
            Msg_AsFem_Exit();
        }
        getline(in,str);linenum+=1;
    }
    return HasType;
}