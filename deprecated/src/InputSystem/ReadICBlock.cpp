#include "InputSystem/InputSystem.h"

bool InputSystem::ReadICBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,ICSystem &icSystem)
{
    // each bc block should looks like:
    //   [ic1]
    //     type=const [random,user1,user2,user3...]
    //     dof=u1
    //     value=1.0 [default is 0.0, so it is not necessary to be given!!!]
    //     block=block_name [i.e. all]
    //   [end]
    // important: now , str already contains [bcs] !!!

    bool HasICBlock=false;
    ICBlock icblock;
    string tempstr,str0;
    vector<double> number;

    bool HasElmt=false;
    bool HasValue=false;
    bool HasBlock=false;
    bool HasDof=false;

    // now str="[bcs]"
    while (linenum<=lastendlinenum)
    {
        getline(in,str);linenum+=1;
        str=RemoveStrSpace(str);
        if(IsCommentLine(str)||str.size()<1) continue;

        HasElmt=false;
        HasValue=false;
        HasBlock=false;
        HasDof=false;

        if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
          (str.find("[end]")==string::npos&&str.find("[END]")==string::npos))
        {
            HasICBlock=false;
            tempstr=RemoveStrSpace(str);
            if(tempstr.size()<3)
            {
                Msg_Input_LineError(linenum);
                Msg_Input_ICSubBlockNameNotFound();
                Msg_AsFem_Exit();
                return false;
            }
            else
            {
                icblock.Clear();
                icblock._ICBlockName=tempstr.substr(1,tempstr.size()-1-1);
                HasICBlock=true;
            }
            while(str.find("[end]")==string::npos&&str.find("[END]")==string::npos)
            {
                getline(in,str);linenum+=1;
                str0=str;
                str=RemoveStrSpace(str0);
                if(IsCommentLine(str)||str.size()<1) continue;
                if(str.find("type=")!=string::npos)
                {
                    string substr=str.substr(str.find('=')+1);
                    if(substr.find("const")!=string::npos&&substr.length()==5)
                    {
                        icblock._ICElmtName="const";
                        icblock._ICTypeID=ICType::CONST;
                        HasElmt=true;
                    }
                    else if(substr.find("random")!=string::npos&&substr.length()==6)
                    {
                        icblock._ICElmtName="random";
                        icblock._ICTypeID=ICType::RAND;
                        HasElmt=true;
                    }
                    else if(substr.find("circle")!=string::npos&&substr.length()==6)
                    {
                        icblock._ICElmtName="circle";
                        icblock._ICTypeID=ICType::CIRCLE;
                        HasElmt=true;
                    }
                    else if(substr.find("user")!=string::npos)
                    {
                        number=SplitStrNum(str);
                        if(number.size()<1)
                        {
                            Msg_Input_LineError(linenum);
                            Msg_Input_ICUserNumNotFound();
                            Msg_AsFem_Exit();
                            return false;
                        }
                        
                        if(int(number[0])<1||int(number[0])>15)
                        {
                            Msg_Input_LineError(linenum);
                            Msg_Input_ICUserNumNotValid();
                            Msg_AsFem_Exit();
                            return false;
                        }
                        else
                        {
                            icblock._ICElmtName=str.substr(str.find_first_of('=')+1);
                            switch(int(number[0])){
                                case 1:
                                    icblock._ICTypeID=ICType::USER1;
                                    break;
                                case 2:
                                    icblock._ICTypeID=ICType::USER2;
                                    break;
                                case 3:
                                    icblock._ICTypeID=ICType::USER3;
                                    break;
                                case 4:
                                    icblock._ICTypeID=ICType::USER4;
                                    break;
                                case 5:
                                    icblock._ICTypeID=ICType::USER5;
                                    break;
                            }
                            HasElmt=true;
                        }
                    }
                    else
                    {
                        Msg_Input_ICBlockUnsupportedType();
                        printf("***        type= %-20s is invalid !!!       ***\n",substr.c_str());
                        Msg_AsFem_Exit();
                    }
                }
                else if(str.find("dof=")!=string::npos)
                {
                    if(str.size()<5)
                    {
                        Msg_Input_LineError(linenum);
                        Msg_Input_DofNameNotFound();
                        Msg_AsFem_Exit();
                    }
                    else
                    {
                        icblock._DofName=str.substr(4);
                        HasDof=true;
                    }
                }
                else if(str.find("block=")!=string::npos)
                {
                    if(str.size()<7)
                    {
                        Msg_Input_LineError(linenum);
                        Msg_Input_ICBlockNameNotFound();
                        Msg_AsFem_Exit();
                        return false;
                    }
                    else
                    {
                        icblock._BlockName=str.substr(str.find_first_of("=")+1);
                        HasBlock=true;
                    }
                }
                else if(str.find("values=")!=string::npos)
                {
                    number=SplitStrNum(str0);
                    if(number.size()<1)
                    {
                        Msg_Input_LineError(linenum);
                        Msg_Input_ICValueNotFound();
                        Msg_AsFem_Exit();
                        return false;
                    }
                    else
                    {
                        icblock._value=number;
                        HasValue=true;
                    }
                }
            }
            if(HasDof&&HasElmt&&HasICBlock)
            {
                HasICBlock=true;
                if(!HasValue) 
                {
                    icblock._value.clear();
                    icblock._value.push_back(0.0);
                    icblock._value.push_back(1.0);
                    Msg_Input_ICBlockNoValueWarning(icblock._ICBlockName);
                }
                if(!HasBlock){
                    icblock._BlockName="alldomain";
                    HasBlock=true;
                }
                icSystem.AddICBlock(icblock);
            }
            else
            {
                Msg_Input_ICBlockInfoNotComplete();
                printf("***   information is missing in [%-20s]     ***\n",icblock._ICBlockName.c_str());
                if(!HasElmt){
                    Msg_Input_ICBlockNoElmt();
                }

                if(!HasDof){
                    Msg_Input_ICBlockNoDof();
                }
               
                if(!HasBlock){
                    Msg_Input_ICBlockNoBlock();
                }
                return false;
            }
        }
    }

    return true;
}