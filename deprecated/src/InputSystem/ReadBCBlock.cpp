#include "InputSystem/InputSystem.h"


bool InputSystem::ReadBCBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,BCSystem &bcSystem)
{
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
    string tempstr;
    vector<double> number;

    bool HasElmt=false;
    bool HasValue=false;
    bool HasBoundary=false;
    bool HasDof=false;

    bcblock._IsTimeDependent=false;

    // now str="[bcs]"
    while (linenum<=lastendlinenum){
        getline(in,str);linenum+=1;
        str=RemoveStrSpace(str);
        if(IsCommentLine(str)||str.size()<1) continue;

        HasElmt=false;
        HasValue=false;
        HasBoundary=false;
        HasDof=false;

        if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
          (str.find("[end]")==string::npos&&str.find("[END]")==string::npos)){
            HasBCBlock=false;
            tempstr=RemoveStrSpace(str);
            if(tempstr.size()<3){
                Msg_Input_LineError(linenum);
                Msg_Input_BCSubBlockNameNotFound();
                Msg_AsFem_Exit();
                return false;
            }
            else
            {
                bcblock.Clear();
                bcblock._BCBlockName=tempstr.substr(1,tempstr.size()-1-1);
                HasBCBlock=true;
            }
            while(str.find("[end]")==string::npos&&str.find("[END]")==string::npos)
            {
                getline(in,str);linenum+=1;
                str=RemoveStrSpace(str);
                if(IsCommentLine(str)||str.size()<1) continue;
                if(str.find("type=")!=string::npos){
                    string substr=str.substr(str.find_first_of('=')+1);
                    if(substr.find("dirichlet")!=string::npos && substr.length()==9){
                        bcblock._BCElmtName="dirichlet";
                        HasElmt=true;
                    }
                    else if(substr.find("neumann")!=string::npos && substr.length()==7){
                        bcblock._BCElmtName="neumann";
                        HasElmt=true;
                    }
                    else if(substr.find("preset")!=string::npos && substr.length()==6){
                        bcblock._BCElmtName="preset";
                        HasElmt=true;
                    }
                    else if(substr.find("user")!=string::npos){
                        number=SplitStrNum(str);
                        if(number.size()<1)
                        {
                            Msg_Input_LineError(linenum);
                            Msg_Input_BCUserNumNotFound();
                            Msg_AsFem_Exit();
                            return false;
                        }
                        
                        if(int(number[0])<1||int(number[0])>15)
                        {
                            Msg_Input_LineError(linenum);
                            Msg_Input_BCUserNumNotValid();
                            Msg_AsFem_Exit();
                            return false;
                        }
                        else
                        {
                            bcblock._BCElmtName=str.substr(str.find_first_of('=')+1);
                            HasElmt=true;
                        }
                    }
                    else{
                        Msg_Input_BCBlockUnsupportedType();
                        printf("***        type= %-20s is invalid !!!       ***\n",substr.c_str());
                        Msg_AsFem_Exit();
                    }
                }
                else if(str.find("dof=")!=string::npos){
                    if(!HasElmt){
                        cout<<"*** Error: type= is not given yet in [bcs] block !!!       ***"<<endl;
                        cout<<"***        dof= should be given after type=      !!!       ***"<<endl;
                        Msg_AsFem_Exit();
                    }
                    if(str.size()<5){
                        Msg_Input_LineError(linenum);
                        Msg_Input_DofNameNotFound();
                        Msg_AsFem_Exit();
                    }
                    else{
                        bcblock._DofName=str.substr(4);
                        HasDof=true;
                    }
                }
                else if(str.find("boundary=")!=string::npos){
                    if(!HasElmt){
                        cout<<"*** Error: type= is not given yet in [bcs] block   !!!     ***"<<endl;
                        cout<<"***        boundary= should be given after type=  !!!      ***"<<endl;
                        Msg_AsFem_Exit();
                    }
                    if(str.size()<10){
                        Msg_Input_LineError(linenum);
                        Msg_Input_BCBoundaryNameNotFound();
                        Msg_AsFem_Exit();
                        return false;
                    }
                    else{
                        bcblock._BoundaryName=str.substr(str.find_first_of("=")+1);
                        HasBoundary=true;
                    }
                }
                else if(str.find("value=")!=string::npos){
                    if(!HasElmt){
                        cout<<"*** Error: type= is not given yet in [bcs] block   !!!     ***"<<endl;
                        cout<<"***        value= should be given after type=      !!!     ***"<<endl;
                        Msg_AsFem_Exit();
                    }
                    number=SplitStrNum(str);
                    if(number.size()<1){
                        Msg_Input_LineError(linenum);
                        Msg_Input_BCValueNotFound();
                        Msg_AsFem_Exit();
                        return false;
                    }
                    else{
                        bcblock._value=number[0];
                        HasValue=true;
                        int i;
                        i=str.find_first_of("=");
                        string substr=str.substr(i+1);
                        // cout<<"bcblock="<<bcblock._BCBlockName<<", str="<<str<<", substr="<<substr<<endl;
                        if(IsValidExpression(substr)){
                            bcblock._IsTimeDependent=true;
                        }
                        else{
                            bcblock._IsTimeDependent=false;
                        }
                    }
                }
            }
            if(HasDof&&HasBoundary&&HasElmt&&HasBCBlock){
                HasBCBlock=true;
                if(!HasValue) bcblock._value=0.0;
                bcSystem.AddBCBlock(bcblock);
            }
            else
            {
                Msg_Input_BCBlockInfoNotComplete();
                printf("***   information is missing in [%-20s]     ***\n",bcblock._BCBlockName.c_str());
                if(!HasElmt){
                    Msg_Input_BCBlockNoElmt();
                }

                if(!HasDof){
                    Msg_Input_BCBlockNoDof();
                }
               
                if(!HasBoundary){
                    Msg_Input_BCBlockNoBoundary();
                }
                return false;
            }
        }
    }

    return true;
}