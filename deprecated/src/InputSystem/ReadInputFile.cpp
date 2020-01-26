#include "InputSystem/InputSystem.h"

bool InputSystem::ReadInputFile(Mesh &mesh,
                                DofHandler &dofHandler,
                                BCSystem &bcSystem,
                                ICSystem &icSystem,
                                ElmtSystem &elmtSystem,
                                MaterialSystem &mateSystem,
                                QPBlockInfo &qpBlockInfo,
                                LinearSolverBlockInfo &linearSolverBlockInfo,
                                NonLinearSolverBlockInfo &nonlinearSolverBlockInfo,
                                SolutionSystem &solutionSystem,
                                FECtrlInfo &feCtrlInfo){
    ifstream in;
    string str;
    int linenum=0;

    bool HasMeshBlock=false;
    bool HasDofsBlock=false;
    bool HasBCBlock=false;
    bool HasICBlock=false;
    bool HasElmtBlock=false;
    bool HasMateBlock=false;
    bool HasJobBlock=false;
    bool HasQpBlock=false;
    bool HasLinearSolverBlock=false;
    bool HasNonLinearSolverBlock=false;
    bool HasProjectionBlock=false;

    if(_HasInputFileName)
    {
        in.open(_InputFileName.c_str(),ios::in);
        while(!in.is_open())
        {
            Msg_InputFileName_Invalid(_InputFileName);
            cout<<"*** Please enter the input file name:";
            cin>>_InputFileName;
        }
    }
    else
    {
        cout<<"*** Please enter the input file name:";
        cin>>_InputFileName;
        in.open(_InputFileName.c_str(),ios::in);
        while(!in.is_open())
        {
            Msg_InputFileName_Invalid(_InputFileName);
            cout<<"*** Please enter the input file name:";
            cin>>_InputFileName;
        }
        _HasInputFileName=true;
    }

    linenum=0;

    HasMeshBlock=false;
    HasDofsBlock=false;
    HasBCBlock=false;
    HasElmtBlock=false;
    HasMateBlock=false;
    HasJobBlock=false;
    HasQpBlock=false;
    HasLinearSolverBlock=false;
    HasNonLinearSolverBlock=false;
    HasProjectionBlock=false;


    while(!in.eof()){
        getline(in,str);linenum+=1;
        str=RemoveStrSpace(str);
        if(IsCommentLine(str)||str.size()<1) continue;
        if(str.find("[Mesh]")!=string::npos||
           str.find("[mesh]")!=string::npos||
           str.find("[MESH]")!=string::npos){
            if(!IsBracketMatch(in,linenum)){
                cout<<"*** Error: [Mesh]/[end] bracket pair not match !!!         ***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
            if(ReadMeshBlock(in,str,linenum,mesh)){
                HasMeshBlock=true;
            }
            else{
                HasMeshBlock=false;
            }
        }
        else if(str.find("[dofs]")!=string::npos||
                str.find("[DOFS]")!=string::npos){
            if(!IsBracketMatch(in,linenum)){
                cout<<"*** Error: [dofs]/[end] bracket pair not match !!!         ***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
            if(ReadDofsBlock(in,str,linenum,dofHandler)){
                HasDofsBlock=true;
            }
            else{
                HasDofsBlock=false;
            }
        }
        else if(str.find("[bcs]")!=string::npos||
                str.find("[BCS]")!=string::npos){
            if(!HasDofsBlock){
                cout<<"*** Error: [bcs] block require [dofs] block   !!!          ***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum))
            {
                if(ReadBCBlock(in,str,lastendlinenum,linenum,bcSystem)){
                    HasBCBlock=true;
                }
                else{
                    HasBCBlock=false;
                }
            }
            else{
                cout<<"*** Error: [bcs]/[end] bracket pair not match !!!          ***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[ics]")!=string::npos||
                str.find("[ICS]")!=string::npos)
        {
            if(!HasDofsBlock){
                cout<<"*** Error: [ics] block require [dofs] block   !!!          ***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadICBlock(in,str,lastendlinenum,linenum,icSystem)){
                    HasICBlock=true;
                }
                else{
                    HasICBlock=false;
                    Msg_AsFem_Exit();
                }
            }
            else{
                cout<<"*** Error: [ics]/[end] bracket pair not match !!!          ***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[elmts]")!=string::npos||
                str.find("[ELMTS]")!=string::npos||
                str.find("[ELMTs]")!=string::npos){
            if(!HasDofsBlock){
                cout<<"*** Error: [elmts] block require [dofs] block !!!          ***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadElmtBlock(in,str,lastendlinenum,linenum,elmtSystem,dofHandler)){
                    HasElmtBlock=true;
                }
                else{
                    HasElmtBlock=false;
                }
            }
            else{
                cout<<"*** Error: [elmts]/[end] bracket pair not match !!!        ***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[mates]")!=string::npos||
                str.find("[MATES]")!=string::npos||
                str.find("[MATEs]")!=string::npos){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadMateBlock(in,str,lastendlinenum,linenum,mateSystem)){
                    HasMateBlock=true;
                }
                else{
                    HasMateBlock=false;
                }
            }
            else{
                cout<<"*** Error: [mates]/[end] bracket pair not match !!!        ***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[job]")!=string::npos||
                str.find("[JOBs]")!=string::npos||
                str.find("[JOBS]")!=string::npos||
                str.find("[jobS]")!=string::npos){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadJobBlock(in,str,linenum,feCtrlInfo)){
                    HasJobBlock=true;
                }
                else{
                    HasJobBlock=false;
                }
            }
            else{
                cout<<"*** Error: [job]/[end] bracket pair not match !!!          ***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[qpoint]")!=string::npos||
                str.find("[QPOINT]")!=string::npos||
                str.find("[Qpoint]")!=string::npos||
                str.find("[QPoint]")!=string::npos){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadQpBlock(in,str,linenum,qpBlockInfo)){
                    HasQpBlock=true;
                }
                else{
                    HasQpBlock=false;
                }
            }
            else{
                cout<<"*** Error: [qpoint]/[end] bracket pair not match !!!       ***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if((str.find("[linearsolver]")!=string::npos||
                str.find("[LinearSolver]")!=string::npos||
                str.find("[LINEARSOLVER]")!=string::npos)&&str.length()==14){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadLinearSolverBlock(in,str,linenum,linearSolverBlockInfo)){
                    HasLinearSolverBlock=true;
                }
                else{
                    HasLinearSolverBlock=false;
                }
            }
            else{
                cout<<"*** Error: [linearsolver]/[end] bracket pair not match !!! ***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if((str.find("[nonlinearsolver]")!=string::npos||
                str.find("[NonLinearSolver]")!=string::npos||
                str.find("[NONLINEARSOLVER]")!=string::npos)&&str.length()==17){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadNonLinearSolverBlock(in,str,linenum,nonlinearSolverBlockInfo)){
                    HasNonLinearSolverBlock=true;
                }
                else{
                    HasNonLinearSolverBlock=false;
                }
            }
            else{
                cout<<"*** Error: [nonlinearsolver]/[end] bracket pair not match !***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if((str.find("[projection]")!=string::npos||
                str.find("[Projection]")!=string::npos||
                str.find("[PROJECTION]")!=string::npos)&&str.length()==12){
            int lastendlinenum;
            if(IsBracketMatch(in,linenum,lastendlinenum)){
                if(ReadProjectionBlock(in,str,linenum,solutionSystem)){
                    HasProjectionBlock=true;
                }
                else{
                    HasProjectionBlock=false;
                }
            }
            else{
                cout<<"*** Error: [projection]/[end] bracket pair not match !!!   ***"<<endl;
                Msg_AsFem_Exit();
                return false;
            }
        }
        else if(str.find("[]")!=string::npos){
            Msg_Input_LineError(linenum);
            Msg_Input_BlockBracketNotComplete();
            Msg_AsFem_Exit();
        }
    }

    if(!HasMeshBlock){
        Msg_InputFile_NoMeshBlock();   
        return false;
    }
    if(!HasDofsBlock){
        Msg_InputFile_NoDofsBlock();
        return false;
    }

    if(!HasBCBlock){
        Msg_Input_NoBCBlockWarning();
    }

    if(!HasElmtBlock){
        Msg_InputFile_NoElmtBlockFound();
        return false;
    }

    if(!HasMateBlock){
        Msg_InputFile_NoMateBlockFoundWarning();
    }

    if(!HasICBlock){
        Msg_InputFile_NoICBlockFoundWarning();
    }

    if(!HasQpBlock){
        qpBlockInfo.Reset();
    }

    if(!HasLinearSolverBlock){
        linearSolverBlockInfo.Reset();
    }

    if(!HasNonLinearSolverBlock){
        nonlinearSolverBlockInfo.Reset();
    }

    if(HasProjectionBlock){};// to get rid of compiler's complain

    if(!HasJobBlock){
        Msg_InputFile_NoJobBlockFound();
        return false;
    }
    return true;
}