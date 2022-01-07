//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.07.08
//+++ Purpose: Implement the reader for [mesh] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/SingleBlockReader.h"
#include "Mesh/Mesh.h"
#include "Mesh/MeshIO.h"

/**
 * Implement the [mesh] block reader class,
 * here the reader will read the mesh information and related options from [mesh] block.
 * If users have no idea about the usage, then just simply using 'type=helper'
 */
class MeshBlockReader:public SingleBlockReader{
public:
    MeshBlockReader(){
        _InputFileName.clear();
        _MeshFileName.clear();
    }
    /**
     * Implement the mesh information reader function, if something is wrong or some information
     * is missing, then it will return false value
     */
    bool ReadMeshBlock(ifstream &in,string str,int &linenum,Mesh &mesh);

    /**
     * This function allow the inputsystem to set the input file name to meshblockreader class,
     * then the mechblockreader can set up the output mesh file name
     * @param name the name of the input file
     * For example, if name='test.i', then the output mesh file name should be 'test_mesh.vtu'
     */
    void SetMeshBlockReaderInputFileName(string name){_InputFileName=name;} 
    /**
     * Implement the helper function for [mesh] block
     */
    virtual void PrintHelper() override;
private:
    string _InputFileName;/**< this string store the input file name, it should be initiliazed in the InputSystem class */
    string _MeshFileName; /**< this string store the name for the output mesh file, defaul format is vtu */

};
