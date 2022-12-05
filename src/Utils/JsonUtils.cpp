//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.31
//+++ Purpose: Implement a general json file reading/element access
//+++          for material properties and other calculations
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Utils/JsonUtils.h"

JsonUtils::JsonUtils(){}

void JsonUtils::checkValidation(const nlohmann::json &t_json,const string &matename){
    if(!t_json.contains(matename)){
        MessagePrinter::printErrorTxt("can\'t find material property("+matename+") in the given json file, please check your input file");
        MessagePrinter::exitAsFem();
    }
}
double JsonUtils::getValue(const nlohmann::json &t_json,const string &matename){
    if(t_json.contains(matename)){
        if(t_json.at(matename).is_number()||
           t_json.at(matename).is_number_float()||
           t_json.at(matename).is_number_integer()||
           t_json.at(matename).is_number_unsigned()){
            return static_cast<double>(t_json.at(matename));
        }
        else{
            MessagePrinter::printErrorTxt("the value of material property(\'"+matename+"\') is not a valid number, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        MessagePrinter::printErrorTxt("can\'t find material property(\'"+matename+"\') in the given json file, please check your input file");
        MessagePrinter::exitAsFem();
    }
    return 0.0;
}

int JsonUtils::getInteger(const nlohmann::json &t_json,const string &matename){
    if(t_json.contains(matename)){
        if(t_json.at(matename).is_number_integer()||
           t_json.at(matename).is_number_unsigned()){
            return static_cast<int>(t_json.at(matename));
        }
        else{
            MessagePrinter::printErrorTxt("the value of material property(\'"+matename+"\') is not a valid integer, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        MessagePrinter::printErrorTxt("can\'t find material property(\'"+matename+"\') in the given json file, please check your input file");
        MessagePrinter::exitAsFem();
    }
    return -1;
}

string JsonUtils::getString(const nlohmann::json &t_json,const string &matename){
    if(t_json.contains(matename)){
        if(t_json.at(matename).is_string()){
            return t_json.at(matename);
        }
        else{
            MessagePrinter::printErrorTxt("the value of material property(\'"+matename+"\') is not a valid string, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        MessagePrinter::printErrorTxt("can\'t find material property(\'"+matename+"\') in the given json file, please check your input file");
        MessagePrinter::exitAsFem();
    }
    return "";
}

bool JsonUtils::getBoolean(const nlohmann::json &t_json,const string &matename){
    if(t_json.contains(matename)){
        if(t_json.at(matename).is_boolean()){
            return t_json.at(matename);
        }
        else{
            MessagePrinter::printErrorTxt("the value of material property(\'"+matename+"\') is not a valid boolean, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        MessagePrinter::printErrorTxt("can\'t find material property(\'"+matename+"\') in the given json file, please check your input file");
        MessagePrinter::exitAsFem();
    }
    return false;
}

Vector3d JsonUtils::getVector(const nlohmann::json &t_json,const string &matename){
    if(t_json.contains(matename)){
        if(!t_json.at(matename).is_array()){
            MessagePrinter::printErrorTxt("the vector value of material property(\'"+matename+"\') is not a valid vector(3d), please check your input file");
            MessagePrinter::exitAsFem();
        }
        else{
            Vector3d temp(0.0);
            for(int i=0;i<static_cast<int>(t_json.at(matename).size());i++){
                if(t_json.at(matename).at(i).is_number()){
                    temp(i+1)=static_cast<double>(t_json.at(matename).at(i));
                }
                else{
                    MessagePrinter::printErrorTxt(to_string(i+1)+"-th value of material property(\'"+matename+"\') is not a valid number, please check your input file");
                    MessagePrinter::exitAsFem();
                }
            }
            return temp;
        }
    }
    else{
        MessagePrinter::printErrorTxt("can\'t find material property(\'"+matename+"\') in the given json file, please check your input file");
        MessagePrinter::exitAsFem();
    }
    return Vector3d(0);
}

bool JsonUtils::hasValue(const nlohmann::json &t_json,const string &matename){
    if(t_json.contains(matename)){
        return true;
    }
    return false;
}


bool JsonUtils::hasOnlyGivenValues(const nlohmann::json &t_json,const vector<string> &namevec){
    bool isValid=true;
    string jsonname;
    for(auto it=t_json.begin();it!=t_json.end();it++){
        isValid=false;
        jsonname=it.key();
        for(int j=0;j<static_cast<int>(namevec.size());j++){
            if(jsonname==namevec[j]){
                isValid=true;break;
            }
        }
        if(!isValid){
            return isValid;
        }
    }
    return isValid;
}