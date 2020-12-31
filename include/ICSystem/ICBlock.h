//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.11
//+++ Purpose: Define [ics] sub block for our input file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include <vector>


#include "ICSystem/ICType.h"

using namespace std;

class ICBlock{
public:
    ICBlock(){
        _ICBlockName.clear();
        _ICTypeName.clear();
        _ICType=ICType::NULLIC;
        _DofName.clear();
        _DofID=-1;
        _Parameters.clear();
        _DomainNameList.clear();// support multiple boundary name list
    }

    string         _ICBlockName;
    string         _ICTypeName;
    ICType         _ICType;
    string         _DofName;
    int            _DofID;
    vector<double> _Parameters;
    vector<string> _DomainNameList;// support multiple boundary name list

    void Init(){
        _ICBlockName.clear();
        _ICTypeName.clear();
        _ICType=ICType::NULLIC;
        _DofName.clear();
        _DofID=-1;
        _Parameters.clear();
        _DomainNameList.clear();// support multiple boundary name list
    }
};