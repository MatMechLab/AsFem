//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
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
#include "SolutionSystem/SolutionSystem.h"


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
    void SetOutputType(OutputType outputtype);
    //************************************************************
    //*** basic getting functions
    //************************************************************
    inline int GetIntervalNum()const{return _Interval;}
    inline string GetOutputFileName()const{return _OutputFileName;}
    inline string GetPVDFileName()const{return _PVDFileName;}
    //************************************************************
    //*** write out our results to files with different format
    //************************************************************
    void WriteResultToFile(const Mesh &mesh,const DofHandler &dofHandler,const SolutionSystem &solutionSystem);
    void WriteResultToFile(const int &step,const Mesh &mesh,const DofHandler &dofHandler,const SolutionSystem &solutionSystem);

    //*************************************************
    //*** for pvd file
    //*************************************************
    void WritePVDFileHeader();
    void WritePVDFileEnd();
    void WriteResultToPVDFile(const double &timestep,string resultfilename);

    void PrintInfo()const;
    
private:
    void WriteResult2VTU(const Mesh &mesh,const DofHandler &dofHandler,const SolutionSystem &solutionSystem);
    void WriteResult2VTU(const int &step,const Mesh &mesh,const DofHandler &dofHandler,const SolutionSystem &solutionSystem);
    void WriteResult2VTK(const Mesh &mesh,const DofHandler &dofHandler,const SolutionSystem &solutionSystem);
    void WriteResult2CSV(const Mesh &mesh,const DofHandler &dofHandler,const SolutionSystem &solutionSystem);
    //**************************

private:
    int _Interval;
    OutputType _OutputType;
    string _OutputFileName,_InputFileName;
    string _OutputTypeName;
    string _OutputFolderName;
    string _PVDFileName;
    vector<string> _CSVFieldNameList;

private:
    //****************************************
    //*** for PETSc vec
    //****************************************
    Vec _Useq,_ProjSeq,_ProjScalarSeq,_ProjVectorSeq,_ProjRank2Seq,_ProjRank4Seq;
    VecScatter _scatterU,_scatterProj,_scatterProjScalar;
    VecScatter _scatterProjVector,_scatterProjRank2,_scatterProjRank4;
    PetscMPIInt _rank;

private:
    //****************************************
    //*** for file I/O
    //****************************************
    string _VTUFileName;
    string _OutputFilePrefix;

};