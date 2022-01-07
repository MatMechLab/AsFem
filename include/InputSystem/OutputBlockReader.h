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
//+++ Date   : 2021.08.05
//+++ Purpose: Implement the reader for [output] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/SingleBlockReader.h"

#include "OutputSystem/OutputSystem.h"

/**
 * This class implement the reader for [qpoint] block
 */
class OutputBlockReader:public SingleBlockReader{
public:
    OutputBlockReader(){
        _InputFileName.clear();
    }

    /**
     * Set up the input file name within the InputSystem class
     */
    void SetOutputBlockReaderInputFileName(string filename){_InputFileName=filename;}

    /**
     * Print out the basic information for the [dofs] block <br>
     * If you dont know what to do, just using 'type=helper' <br>
     */
    virtual void PrintHelper() override;

    /**
     * Read the [output] block, if everything works fine, it should return True value<br>
     * The block should looks like: <br>
     * <pre>
     * [output]
     *   type=vtu[vtk,csv,txt]
     *   folder=foldername[default is empty]
     *   interval=2
     * [end]
     * </pre>
     * @param in ifstream for the file reading
     * @param str the string which constains '[dofs]' line
     * @param linenum the current line number
     * @param outputSystem OutputSystem class, which is used to set up the basic information for the output
     */
    bool ReadOutputBlock(ifstream &in,string str,int &linenum,OutputSystem &outputSystem);

private:
    string _InputFileName;/**< store the string name of the input file*/

};
