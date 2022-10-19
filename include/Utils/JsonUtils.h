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

#pragma once

#include "Utils/MessagePrinter.h"
#include "nlohmann/json.hpp"
#include "MathUtils/Vector3d.h"


/**
 * This class implement the general parameter access to json file for material
 * properties calculation and other operations
 */
class JsonUtils{
public:
    /**
     * constructor
     */
    JsonUtils();

    /**
     * check wether the parameter is exist
     * @param t_json the input json class
     * @param matename the string name of the material property
     */
    static void checkValidation(const nlohmann::json &t_json,const string &matename);

    /**
     * get the parameter's value via its string name
     * @param t_json the input json class
     * @param matename the string name of the material property
     */
    static double getValue(const nlohmann::json &t_json,const string &matename);

    /**
     * get the parameter's integer value via its string name
     * @param t_json the input json class
     * @param matename the string name of the material property
     */
    static int getInteger(const nlohmann::json &t_json,const string &matename);
    
    /**
     * get the parameter's boolean value via its string name
     * @param t_json the input json class
     * @param matename the string name of the material property
     */
    static bool getBoolean(const nlohmann::json &t_json,const string &matename);

    /**
     * get the parameter's string value via its string name
     * @param t_json the input json class
     * @param matename the string name of the material property
     */
    static string getString(const nlohmann::json &t_json,const string &matename);

    /**
     * get the parameter's vector value via its string name
     * @param t_json the input json class
     * @param matename the string name of the material property
     */
    static Vector3d getVector(const nlohmann::json &t_json,const string &matename);

    /**
     * check wethear the value name is valid
     * @param t_json the input json class
     * @param matename the string name of the material property
     */
    static bool hasValue(const nlohmann::json &t_json,const string &matename);

    /**
     * check wether the json only contains the given variables
     * @param t_json the input json class
     * @param namevec the string name vector of the material property
     */
    static bool hasOnlyGivenValues(const nlohmann::json &t_json,const vector<string> &namevec);

};