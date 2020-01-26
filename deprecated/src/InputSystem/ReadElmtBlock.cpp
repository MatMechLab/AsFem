#include "InputSystem/InputSystem.h"

bool InputSystem::ReadElmtBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,ElmtSystem &elmtSystem,DofHandler &dofHandler)
{
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

    // now the str is [elmts]
    while(linenum<=lastendlinenum)
    {
        getline(in,str);linenum+=1;
        str0=str;
        str=RemoveStrSpace(str);
        if(IsCommentLine(str)||str.size()<1) continue;

        HasElmtBlock=false;
        HasElmtType=false;
        HasDofs=false;
        HasBlock=false;
        HasMate=false;

        if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
          (str.find("[end]")==string::npos&&str.find("[END]")==string::npos))
        {
            HasElmtBlock=false;
            tempstr=RemoveStrSpace(str);
            if(tempstr.size()<3)
            {
                Msg_Input_LineError(linenum);
                Msg_Input_ElmtSubBlockNameNotFound();
                Msg_AsFem_Exit();
                return false;
            }
            else
            {
                elmtBlock.Clean();
                elmtBlock._ElmtBlockName=tempstr.substr(tempstr.find_first_of('[')+1,tempstr.find_first_of(']')-1);
                HasElmtBlock=true;
            }
            while(str.find("[end]")==string::npos&&str.find("[END]")==string::npos)
            {
                getline(in,str);linenum+=1;
                str0=str;
                str=RemoveStrSpace(str);

                if(IsCommentLine(str)||str.size()<1) 
                {
                    //cout<<"str="<<str<<endl;
                    continue;
                }

                if(str.find("type=")!=string::npos)
                {
                    substr=str.substr(str.find_first_of('=')+1);
                    if(substr.find("poisson")!=string::npos && substr.length()==7)
                    {
                        elmtBlock._ElmtTypeName="poisson";
                        elmtBlock._ElmtTypeID=ElmtType::POISSON;
                        HasElmtType=true;
                    }
                    else if(substr.find("diffusion")!=string::npos && substr.length()==9)
                    {
                        elmtBlock._ElmtTypeName="diffusion";
                        elmtBlock._ElmtTypeID=ElmtType::DIFFUSION;
                        HasElmtType=true;
                    }
                    else if(substr.find("wave")!=string::npos && substr.length()==4)
                    {
                        elmtBlock._ElmtTypeName="wave";
                        elmtBlock._ElmtTypeID=ElmtType::WAVE;
                        HasElmtType=true;
                    }
                    else if(substr.find("mechanics")!=string::npos && substr.length()==9)
                    {
                        elmtBlock._ElmtTypeName="mechanics";
                        elmtBlock._ElmtTypeID=ElmtType::MECHANICS;
                        HasElmtType=true;
                    }
                    else if(substr.find("cahnhilliard")!=string::npos && substr.length()==12)
                    {
                        elmtBlock._ElmtTypeName="cahnhilliard";
                        elmtBlock._ElmtTypeID=ElmtType::CAHNHILLIARD;
                        HasElmtType=true;
                    }
                    else if(substr.find("cahnhilliardmechanics")!=string::npos && substr.length()==21)
                    {
                        elmtBlock._ElmtTypeName="cahnhilliardmechanics";
                        elmtBlock._ElmtTypeID=ElmtType::CAHNHILLIARDMECHANICS;
                        HasElmtType=true;
                    }
                    else if(substr.find("thermalmechanics")!=string::npos && substr.length()==16)
                    {
                        elmtBlock._ElmtTypeName="thermalmechanics";
                        elmtBlock._ElmtTypeID=ElmtType::THERMALMECHANICS;
                        HasElmtType=true;
                    }
                    else if(substr.find("linearelasticphasefieldfracture")!=string::npos && substr.length()==31)
                    {
                        elmtBlock._ElmtTypeName="linearelasticphasefieldfracture";
                        elmtBlock._ElmtTypeID=ElmtType::LINEARELASTICPHASEFIELDFRACTURE;
                        HasElmtType=true;
                    }
                    else if(substr.find("user")!=string::npos)
                    {
                        number=SplitStrNum(substr);
                        if(number.size()<1)
                        {
                            Msg_Input_LineError(linenum);
                            Msg_Input_ElmtUserNumNotFound();
                            Msg_AsFem_Exit();
                            return false;
                        }
                        if(int(number[0])<1||int(number[0])>20)
                        {
                            Msg_Input_LineError(linenum);
                            Msg_Input_ElmtUserNumNotValid();
                            Msg_AsFem_Exit();
                            return false;
                        }
                        else
                        {
                            int iuel=int(number[0]);
                            elmtBlock._ElmtTypeName=str.substr(str.find_first_of('=')+1);
                            if(iuel==1){
                                elmtBlock._ElmtTypeID=ElmtType::USER1;
                            }
                            else if(iuel==2){
                                elmtBlock._ElmtTypeID=ElmtType::USER2;
                            }
                            else if(iuel==3){
                                elmtBlock._ElmtTypeID=ElmtType::USER3;
                            }
                            else if(iuel==4){
                                elmtBlock._ElmtTypeID=ElmtType::USER4;
                            }
                            else if(iuel==5){
                                elmtBlock._ElmtTypeID=ElmtType::USER5;
                            }
                            else if(iuel==6){
                                elmtBlock._ElmtTypeID=ElmtType::USER6;
                            }
                            else if(iuel==7){
                                elmtBlock._ElmtTypeID=ElmtType::USER7;
                            }
                            else if(iuel==8){
                                elmtBlock._ElmtTypeID=ElmtType::USER8;
                            }
                            else if(iuel==9){
                                elmtBlock._ElmtTypeID=ElmtType::USER9;
                            }
                            else if(iuel==10){
                                elmtBlock._ElmtTypeID=ElmtType::USER10;
                            }
                            else if(iuel==11){
                                elmtBlock._ElmtTypeID=ElmtType::USER11;
                            }
                            else if(iuel==12){
                                elmtBlock._ElmtTypeID=ElmtType::USER12;
                            }
                            else if(iuel==13){
                                elmtBlock._ElmtTypeID=ElmtType::USER13;
                            }
                            else if(iuel==14){
                                elmtBlock._ElmtTypeID=ElmtType::USER14;
                            }
                            else if(iuel==15){
                                elmtBlock._ElmtTypeID=ElmtType::USER15;
                            }
                            else if(iuel==16){
                                elmtBlock._ElmtTypeID=ElmtType::USER16;
                            }
                            else if(iuel==17){
                                elmtBlock._ElmtTypeID=ElmtType::USER17;
                            }
                            else if(iuel==18){
                                elmtBlock._ElmtTypeID=ElmtType::USER18;
                            }
                            else if(iuel==19){
                                elmtBlock._ElmtTypeID=ElmtType::USER19;
                            }
                            else if(iuel==20){
                                elmtBlock._ElmtTypeID=ElmtType::USER20;
                            }
                            HasElmtType=true;
                        }
                    }
                    else
                    {
                        Msg_Input_ElmtBlockUnsupportedType();
                        printf("***        type= %-20s is invalid !!!       ***\n",substr.c_str());
                        Msg_AsFem_Exit();
                        return false;
                    }
                }
                else if(str.find("dofs=")!=string::npos)
                {
                    if(str.size()<5)
                    {
                        Msg_Input_LineError(linenum);
                        Msg_Input_DofNameNotFound();
                        Msg_AsFem_Exit();
                    }
                    else
                    {
                        substr=str0.substr(str0.find_first_of('=')+1);
                        vector<string> strtmp=SplitStr(substr,' ');
                        if(IsUniqueStrVec(strtmp)){
                            elmtBlock._ElmtDofNameList=SplitStr(substr,' ');
                            elmtBlock._ElmtDofIndexList=dofHandler.GetDofsIndexFromNameVec(elmtBlock._ElmtDofNameList);
                            HasDofs=true;
                        }
                        else{
                            HasDofs=false;
                            Msg_Input_LineError(linenum);
                            Msg_Input_ElmtBlockDuplicateDofsName();
                            Msg_AsFem_Exit();
                        }
                    }
                }
                else if(str.find("mate=")!=string::npos)
                {
                    if(str.size()<5)
                    {
                        Msg_Input_LineError(linenum);
                        Msg_Input_ElmtBlockMateNameNotFound();
                        Msg_AsFem_Exit();
                        return false;
                    }
                    else
                    {
                        elmtBlock._ElmtMateName=str.substr(str.find_first_of("=")+1);
                        HasMate=true;
                    }
                }
                else if(str.find("block=")!=string::npos)
                {
                    if(str.size()<7)
                    {
                        Msg_Input_LineError(linenum);
                        Msg_Input_ElmtBlockBlockNameNotFound();
                        Msg_AsFem_Exit();
                        return false;
                    }
                    else
                    {
                        elmtBlock._ElmtDomainBlockName=str.substr(str.find_first_of("=")+1);
                        HasBlock=true;
                    }
                }
            }
            if(HasElmtBlock&&HasDofs&&HasElmtType)
            {
                if(!HasMate){
                    elmtBlock._ElmtMateName="";
                    elmtBlock._MateTypeID=MaterialType::CONSTPOISSONMATE;
                }
                if(!HasBlock) elmtBlock._ElmtDomainBlockName="alldomain";
                elmtSystem.AddElmtBlock(elmtBlock);
            }
            else
            {
                Msg_Input_ElmtBlockInfoNotComplete();
                printf("***   information is missing in [%-20s]     ***\n",elmtBlock._ElmtBlockName.c_str());
                if(!HasElmtType)
                {
                    Msg_Input_ElmtBlockNoType();
                }

                if(!HasDofs)
                {
                    Msg_Input_ElmtBlockNoDofs();
                }
                return false;
            }
        }
    }

    return true;
}