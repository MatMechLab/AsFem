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
//+++ Date   : 2020.07.12
//+++ Purpose: define the output system for AsFem, where all the 
//+++          results should be written out to the result file
//+++          by this class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

#include "OutputSystem/OutputBlock.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"


using namespace std;


class OutputSystem{
public:
    OutputSystem();

    void Init(string inputfilename);
    void InitFromOutputBlock(OutputBlock &outputblock);
    //************************************************************
    //*** basic settings
    //************************************************************
    inline void SetInputFileName(string inputfilename){_InputFileName=inputfilename;}

    //************************************************************
    //*** basic getting functions
    //************************************************************
    inline int GetIntervalNum()const{return _Interval;}
    inline string GetOutputFileName()const{return _OutputFileName;}

    //************************************************************
    //*** write out our results to files with different format
    //************************************************************
    void WriteResultToFile(const Mesh &mesh,const DofHandler &dofHandler,const Vec &U);
    void WriteResultToFile(const Mesh &mesh,const DofHandler &dofHandler,const Vec &U,const int &nProj,const Vec &Proj);

private:
    void WriteResult2VTU(const Mesh &mesh,const DofHandler &dofHandler,const Vec &U);
    void WriteResult2VTK(const Mesh &mesh,const DofHandler &dofHandler,const Vec &U);
    void WriteResult2CSV(const Mesh &mesh,const DofHandler &dofHandler,const Vec &U);
    //**************************
    void WriteResult2VTU(const Mesh &mesh,const DofHandler &dofHandler,const Vec &U,const int &nProj,const Vec &Proj);
    void WriteResult2VTK(const Mesh &mesh,const DofHandler &dofHandler,const Vec &U,const int &nProj,const Vec &Proj);
    void WriteResult2CSV(const Mesh &mesh,const DofHandler &dofHandler,const Vec &U,const int &nProj,const Vec &Proj);


private:
    int _Interval;
    OutputType _OutputType;
    string _OutputFileName,_InputFileName;
    string _OutputTypeName;
    string _OutputFolderName;
    vector<string> _CSVFieldNameList;

private:
    //****************************************
    //*** for PETSc vec
    //****************************************
    Vec _Useq,_ProjSeq;
    VecScatter _scatterU,_scatterProj;
    PetscMPIInt _rank;

private:
    //****************************************
    //*** for file I/O
    //****************************************
    string _VTUFileName;
    string _OutputFilePrefix;

};