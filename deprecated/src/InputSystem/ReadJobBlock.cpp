#include "InputSystem/InputSystem.h"

bool InputSystem::ReadJobBlock(ifstream &in,string str,int &linenum,
                               FECtrlInfo &feCtrlInfo){
    bool HasType=false;
    bool HasEndTime=false;
    bool HasDt=false;
    vector<double> numbers;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;
    numbers.clear();
    feCtrlInfo._IsPrintMesh=false;
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            continue;
        }
        str=RemoveStrSpace(str);
        
        if(str.find("type=")!=string::npos||
           str.find("TYPE=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if(substr.find("static")!=string::npos||
               substr.find("STATIC")!=string::npos){
                   feCtrlInfo._JobType=FEJobType::STATIC;
                   HasType=true;
            }
            else if(substr.find("transient")!=string::npos||
               substr.find("TRANSIENT")!=string::npos){
                   feCtrlInfo._JobType=FEJobType::TRANSIENT;
                   HasType=true;
            }
            else{
                HasType=false;
                Msg_Input_LineError(linenum);
                Msg_Input_JobBlock_InvalidJobType();
                Msg_AsFem_Exit();
            }
        }
        else if(str.find("debug=")!=string::npos||
                str.find("DEBUG=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if(substr.find("true")!=string::npos||
               substr.find("TRUE")!=string::npos){
                feCtrlInfo._IsDebugOn=true;
                feCtrlInfo._IsDebugDep=false;
            }
            else if(substr.find("false")!=string::npos||
                    substr.find("FALSE")!=string::npos){
                feCtrlInfo._IsDebugOn=false;
                feCtrlInfo._IsDebugDep=false;
            }
            else if(substr.find("dep")!=string::npos||
                    substr.find("DEP")!=string::npos){
                feCtrlInfo._IsDebugOn=true;
                feCtrlInfo._IsDebugDep=true;
            }
            else{
                Msg_Input_LineError(linenum);
                Msg_Input_JobBlock_InvalidDebugOption();
                Msg_AsFem_Exit();
            }
        }
        else if(str.compare(0,11,"timemethod=")==0||str.compare(0,11,"TIMEMETHOD=")==0){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if(substr.find("backwardeuler")!=string::npos||
               substr.find("BACKWARDEULER")!=string::npos){
                   feCtrlInfo._TimeSteppingMethod=TimeSteppingType::BACKWARDEULER;
            }
            else if(substr.find("be")!=string::npos||
                    substr.find("BE")!=string::npos){
                feCtrlInfo._TimeSteppingMethod=TimeSteppingType::BACKWARDEULER;
            }
            else if(substr.find("cranknicolson")!=string::npos||
                    substr.find("CRANKNICOLSON")!=string::npos){
                feCtrlInfo._TimeSteppingMethod=TimeSteppingType::CRANKNICOLSON;
            }
            else if(substr.find("cn")!=string::npos||
                    substr.find("CN")!=string::npos){
                feCtrlInfo._TimeSteppingMethod=TimeSteppingType::CRANKNICOLSON;
            }
            else if(substr.find("bdf2")!=string::npos||
                    substr.find("BDF2")!=string::npos){
                feCtrlInfo._TimeSteppingMethod=TimeSteppingType::BDF2;
            }
            else if(substr.find("bdf3")!=string::npos||
                    substr.find("BDF3")!=string::npos){
                feCtrlInfo._TimeSteppingMethod=TimeSteppingType::BDF3;
            }
            else if(substr.find("bdf4")!=string::npos||
                    substr.find("BDF4")!=string::npos){
                feCtrlInfo._TimeSteppingMethod=TimeSteppingType::BDF4;
            }
            else if(substr.find("bdf5")!=string::npos||
                    substr.find("BDF5")!=string::npos){
                feCtrlInfo._TimeSteppingMethod=TimeSteppingType::BDF5;
            }
            else{
                Msg_Input_LineError(linenum);
                Msg_Input_JobBlock_InvalidTimeSteppingOption();
                Msg_AsFem_Exit();
            }
        }
        else if(str.find("endtime=")!=string::npos||
                str.find("ENDTIME=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_JobBlock_NoTimeFound();
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<1.0e-10||numbers[0]>1.0e20){
                    Msg_Input_LineError(linenum);
                    Msg_Input_JobBlock_InvalidTimeValue();
                    Msg_AsFem_Exit();
                }
                feCtrlInfo._EndTime=numbers[0];
                HasEndTime=true;
            }
        }
        else if(str.find("dt=")!=string::npos||
                str.find("DT=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_JobBlock_NoDtFound();
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<1.0e-10){
                    Msg_Input_LineError(linenum);
                    Msg_Input_JobBlock_InvalidDtValue();
                    Msg_AsFem_Exit();
                }
                feCtrlInfo._dt0=numbers[0];
                feCtrlInfo._dt=feCtrlInfo._dt0;
                HasDt=true;
            }
        }
        else if(str.find("interval=")!=string::npos||
                str.find("INTERVAL=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_JobBlock_NoIntervalFound();
                Msg_AsFem_Exit();
            }
            else{
                if(int(numbers[0])<1){
                    Msg_Input_LineError(linenum);
                    Msg_Input_JobBlock_InvalidIntervalValue();
                    Msg_AsFem_Exit();
                }
                feCtrlInfo._Interval=int(numbers[0]);
                HasDt=true;
            }
        }
        else if(str.find("projection=")!=string::npos||
                str.find("PROJECTION=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if(substr.find("true")!=string::npos||
               substr.find("TRUE")!=string::npos){
                feCtrlInfo._IsProjOn=true;
            }
            else if(substr.find("false")!=string::npos||
                    substr.find("FALSE")!=string::npos){
                feCtrlInfo._IsProjOn=false;
            }
            else{
                Msg_Input_LineError(linenum);
                Msg_Input_JobBlock_InvalidProjectionOption();
                Msg_AsFem_Exit();
            }
        }
        else if(str.find("printmesh=")!=string::npos||
                str.find("Printmesh=")!=string::npos||
                str.find("PrintMesh=")!=string::npos||
                str.find("PrintMesh=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if(substr.find("true")!=string::npos||
               substr.find("TRUE")!=string::npos){
                feCtrlInfo._IsPrintMesh=true;
            }
            else if(substr.find("false")!=string::npos||
                    substr.find("FALSE")!=string::npos){
                feCtrlInfo._IsPrintMesh=false;
            }
            else{
                Msg_Input_LineError(linenum);
                Msg_Input_JobBlock_InvalidProjectionOption();
                Msg_AsFem_Exit();
            }
        }
        else if(str.find("adaptive=")!=string::npos||
                str.find("ADAPTIVE=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if(substr.find("true")!=string::npos||
               substr.find("TRUE")!=string::npos){
                feCtrlInfo._IsAdaptive=true;
            }
            else if(substr.find("false")!=string::npos||
                    substr.find("FALSE")!=string::npos){
                feCtrlInfo._IsAdaptive=false;
            }
            else{
                Msg_Input_LineError(linenum);
                Msg_Input_JobBlock_InvalidAdaptiveOption();
                Msg_AsFem_Exit();
            }
        }
        else if(str.find("dtmax=")!=string::npos||
                str.find("DTMAX=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_JobBlock_NoIntervalFound();
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<1.0e-12){
                    Msg_Input_LineError(linenum);
                    Msg_Input_JobBlock_InvalidIntervalValue();
                    Msg_AsFem_Exit();
                }
                feCtrlInfo._DtMax=numbers[0];
                HasDt=true;
            }
        }
        else if(str.find("dtmin=")!=string::npos||
                str.find("DTMIN=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_JobBlock_NoIntervalFound();
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<1.0e-12){
                    Msg_Input_LineError(linenum);
                    Msg_Input_JobBlock_InvalidIntervalValue();
                    Msg_AsFem_Exit();
                }
                feCtrlInfo._DtMin=numbers[0];
                HasDt=true;
            }
        }
        else if(str.find("optiters=")!=string::npos||
                str.find("Optiters=")!=string::npos||
                str.find("OptIters=")!=string::npos||
                str.find("OPTITERS=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                Msg_Input_JobBlock_NoDtFound();
                Msg_AsFem_Exit();
            }
            else{
                if(int(numbers[0])<1){
                    Msg_Input_LineError(linenum);
                    cout<<"*** Error: invalid optimal iters in [job] block !!!  ***"<<endl;
                    Msg_AsFem_Exit();
                }
                feCtrlInfo._OptIters=int(numbers[0]);
            }
        }
        else if(str.find("[]")!=string::npos){
            Msg_Input_LineError(linenum);
            Msg_Input_BlockBracketNotComplete();
            Msg_AsFem_Exit();
        }
        else{
            Msg_Input_LineError(linenum);
            cout<<"*** Error: unknown option in [job] block        !!!        ***"<<endl;
            Msg_AsFem_Exit();
        }
        getline(in,str);linenum+=1;
    }
    if(feCtrlInfo._JobType==FEJobType::TRANSIENT){
        if(!HasEndTime&&!HasDt){
            Msg_Input_JobBblock_TransientInfoNotComplete();
            Msg_AsFem_Exit();
        }
    }

    return HasType;
}