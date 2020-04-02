//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_OUTPUTSYSTEM_H
#define ASFEM_OUTPUTSYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>

//#include <filesystem>

#include <cstdlib>

#include "petsc.h"

//*********************************
//*** AsFem's own header file
//*********************************
#include "MessagePrinter/MessagePrinter.h"
#include "Utils/StringUtils.h"
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "OutputBlock.h"

using namespace std;



class OutputSystem{
public:
    OutputSystem();
    inline string GetLogFileName() const{return _LogFileName;}
    inline string GetVTUFileName() const{return _VTUFileName;}
    void InitOutputStream();
    void SetInputFileName(string inputname) {_InputFileName=inputname;_IsInit=false;}

    void WriteResultToVTU(Mesh &mesh,DofHandler &dofHandler,const Vec &U);
    void WriteResultToVTU(Mesh &mesh,DofHandler &dofHandler,const Vec &U,const int &nproj,const vector<string> &namelist,const Vec &Proj);
    
    // for time step output
    void WriteResultToVTU(const int &step,Mesh &mesh,DofHandler &dofHandler,const Vec &U);
    void WriteResultToVTU(const int &step,Mesh &mesh,DofHandler &dofHandler,const Vec &U,const int &nproj,const vector<string> &namelist,const Vec &Proj);


    OutputBlock _OutputBlock;
    void PrintOutputSystem() const{
        _OutputBlock.PrintOutputBlock();
    }
private:
    bool _IsInit=false;
    bool _IsLogOn=false;
    string _InputFileName;
    string _OutputFilePrefix;
    string _LogFileName;
    // ofstream _VTUFile,_LogFile;
    string _VTUFileName;
    

    //************************************
    //*** For petsc related variables
    PetscMPIInt _rank;
    VecScatter _scatter,_scatterproj;
    Vec _Useq,_PROJseq;

};



#endif // ASFEM_OUTPUTSYSTEM_H