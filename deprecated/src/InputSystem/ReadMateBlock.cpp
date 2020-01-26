#include "InputSystem/InputSystem.h"

bool InputSystem::ReadMateBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,MaterialSystem &mateSystem)
{
    // each bc block should looks like:
    //   [mate1]
    //     type=elastic
    //     params=1.0 2.0
    //   [end]
    // important: now , str already contains [mate] !!!

    bool HasMateBlock=false;
    MaterialBlock mateBlock;
    string tempstr,str0;
    vector<double> number;

    bool HasElmt=false;
    bool HasValue=false;
    bool HasBlock=false;

    // now str="[mate]"
    while (linenum<=lastendlinenum)
    {
        getline(in,str);linenum+=1;
        str=RemoveStrSpace(str);
        if(IsCommentLine(str)||str.size()<1) continue;

        HasElmt=false;
        HasValue=false;
        HasBlock=false;

        if((str.find('[')!=string::npos&&str.find(']')!=string::npos)&&
          (str.find("[end]")==string::npos&&str.find("[END]")==string::npos))
        {
            HasBlock=false;
            tempstr=RemoveStrSpace(str);
            if(tempstr.size()<3)
            {
                Msg_Input_LineError(linenum);
                Msg_Input_MateSubBlockNameNotFound();
                return false;
            }
            else
            {
                mateBlock.Clean();
                mateBlock._MateBlockName=tempstr.substr(1,tempstr.size()-1-1);
                HasBlock=true;
            }
            while(str.find("[end]")==string::npos&&str.find("[END]")==string::npos)
            {
                getline(in,str);linenum+=1;
                str0=str;
                str=RemoveStrSpace(str0);
                if(IsCommentLine(str)||str.size()<1) continue;
                if(str.find("type=")!=string::npos)
                {
                    string substr=str.substr(str.find_first_of('=')+1);
                    if(substr.find("constpoisson")!=string::npos&&substr.length()==12)
                    {
                        mateBlock._MateTypeName="constpoisson";
                        mateBlock._MateTypeID=MaterialType::CONSTPOISSONMATE;
                        HasElmt=true;
                    }
                    else if(substr.find("nonlinearpoisson")!=string::npos&&substr.length()==16)
                    {
                        mateBlock._MateTypeName="nonlinearpoisson";
                        mateBlock._MateTypeID=MaterialType::NONLINEARPOISSONMATE;
                        HasElmt=true;
                    }
                    else if(substr.find("constdiffusion")!=string::npos&&substr.length()==14)
                    {
                        mateBlock._MateTypeName="constdiffusion";
                        mateBlock._MateTypeID=MaterialType::CONSTDIFFUSIONMATE;
                        HasElmt=true;
                    }
                    else if(substr.find("nonlineardiffusion")!=string::npos&&substr.length()==18)
                    {
                        mateBlock._MateTypeName="nonlineardiffusion";
                        mateBlock._MateTypeID=MaterialType::NONLINEARDIFFUSIONMATE;
                        HasElmt=true;
                    }
                    else if(substr.find("cahnhilliard")!=string::npos&&substr.length()==12){
                        mateBlock._MateTypeName="cahnhilliard";
                        mateBlock._MateTypeID=MaterialType::CAHNHILLIARDMATE;
                        HasElmt=true;
                    }
                    else if(substr.find("linearelastic")!=string::npos&&substr.length()==13){
                        mateBlock._MateTypeName="linearelastic";
                        mateBlock._MateTypeID=MaterialType::LINEARELASTICMATE;
                        HasElmt=true;
                    }
                    else if(substr.find("finitestrain")!=string::npos&&substr.length()==12){
                        mateBlock._MateTypeName="finitestrain";
                        mateBlock._MateTypeID=MaterialType::FINITESTRAINMATE;
                        HasElmt=true;
                    }
                    else if(substr.find("saintvenant")!=string::npos&&substr.length()==11){
                        mateBlock._MateTypeName="saintvenant";
                        mateBlock._MateTypeID=MaterialType::STVENANTMATE;
                        HasElmt=true;
                    }
                    else if(substr.find("neohookean")!=string::npos&&substr.length()==10){
                        mateBlock._MateTypeName="neohookean";
                        mateBlock._MateTypeID=MaterialType::NEOHOOKEANMATE;
                        HasElmt=true;
                    }
                    else if(substr.find("generalwave")!=string::npos&&substr.length()==11){
                        mateBlock._MateTypeName="generalwave";
                        mateBlock._MateTypeID=MaterialType::GENERALWAVEMATE;
                        HasElmt=true;
                    }
                    else if(substr.find("thermelastic")!=string::npos&&substr.length()==12){
                        mateBlock._MateTypeName="thermelastic";
                        mateBlock._MateTypeID=MaterialType::ELASTICTHERMALMEHCANICSMATE;
                        HasElmt=true;
                    }
                    else if(substr.find("lpffracture")!=string::npos&&substr.length()==11){
                        mateBlock._MateTypeName="lpffracture";
                        mateBlock._MateTypeID=MaterialType::LINEARELASTICPHASEFIELDFRACTUREMATE;
                        HasElmt=true;
                    }
                    else if(substr.find("modifyfrac")!=string::npos&&substr.length()==10){
                        mateBlock._MateTypeName="modifyfrac";
                        mateBlock._MateTypeID=MaterialType::MODIFYELASTICFRACTUREMATE;
                        HasElmt=true;
                    }
                    else if(substr.find("user")!=string::npos){
                        number=SplitStrNum(substr);
                        if(number.size()<1){
                            Msg_Input_LineError(linenum);
                            Msg_Input_MateUserNumNotFound();
                            Msg_AsFem_Exit();
                            return false;
                        }
                        
                        if(int(number[0])<1||int(number[0])>15){
                            Msg_Input_LineError(linenum);
                            Msg_Input_MateUserNumNotValid();
                            Msg_AsFem_Exit();
                            return false;
                        }
                        else
                        {
                            mateBlock._MateTypeName=str.substr(str.find_first_of('=')+1);
                            switch(int(number[0])){
                                case 1:
                                  mateBlock._MateTypeID=MaterialType::USER1MATE;
                                  break;
                                case 2:
                                  mateBlock._MateTypeID=MaterialType::USER2MATE;
                                  break;
                                case 3:
                                  mateBlock._MateTypeID=MaterialType::USER3MATE;
                                  break;
                                case 4:
                                  mateBlock._MateTypeID=MaterialType::USER4MATE;
                                  break;
                                case 5:
                                  mateBlock._MateTypeID=MaterialType::USER5MATE;
                                  break;
                                case 6:
                                  mateBlock._MateTypeID=MaterialType::USER6MATE;
                                  break;
                                case 7:
                                  mateBlock._MateTypeID=MaterialType::USER7MATE;
                                  break;
                                case 8:
                                  mateBlock._MateTypeID=MaterialType::USER8MATE;
                                  break;
                                case 9:
                                  mateBlock._MateTypeID=MaterialType::USER9MATE;
                                  break;
                                case 10:
                                  mateBlock._MateTypeID=MaterialType::USER10MATE;
                                  break;
                                default:
                                  mateBlock._MateTypeID=MaterialType::NULLMATE;
                                  break;
                            }
                            HasElmt=true;
                        }
                    }
                    else
                    {
                        Msg_Input_MateBlockUnsupportType();
                        printf("***        type= %-20s is invalid !!!       ***\n",substr.c_str());
                        Msg_AsFem_Exit();
                    }
                }
                else if(str0.find("params=")!=string::npos)
                {
                    if(!HasElmt){
                        cout<<"*** Error: params must be given after type in [mates]!     ***"<<endl;
                        Msg_AsFem_Exit();
                    }
                    number=SplitStrNum(str0);
                    if(number.size()<1)
                    {
                        Msg_Input_LineError(linenum);
                        Msg_Input_MateParamsNotFound();
                        return false;
                    }
                    else
                    {
                        mateBlock._MateParams.clear();
                        mateBlock._MateParams=number;
                        HasValue=true;
                    }
                }
            }
            if(HasBlock&&HasElmt)
            {
                HasBlock=true;
                if(!HasValue) 
                {
                    mateBlock._MateParams.clear();
                    mateBlock._MateParams.push_back(1.0);
                    Msg_Input_MateBlockNoParamsWarning(mateBlock._MateBlockName);
                }
                mateSystem.AddMaterialBlock(mateBlock);
                HasMateBlock=true;
            }
            else
            {
                Msg_Input_MateBlockInfoNotComplete();
                printf("***   information is missing in [%-20s]     ***\n",mateBlock._MateBlockName.c_str());
                if(!HasElmt)
                {
                    Msg_Input_MateBlockNoElmt();
                }
                if(!HasValue)
                {
                    Msg_Input_MateBlockNoParams();
                }
                return false;
            }
        }
    }

    return HasMateBlock;
}