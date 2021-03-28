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
//+++ Date   : 2021.02.21
//+++ Purpose: implement the post-process system in AsFem
//+++          this class can do the volume,surface,nodal integration
//+++          or other types post process, and write out the result
//+++          to the csv file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "Postprocess/PostprocessBlock.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "FE/FE.h"
#include "SolutionSystem/SolutionSystem.h"

#include "petsc.h"

using namespace std;

class Postprocess{
public:
    Postprocess();

    void SetInputFileName(string inputfilename);
    void SetOutputInterval(const int &interval){_OutputInterval=interval;}
    inline int GetOutputIntervalNum()const{return _OutputInterval;}
    void AddPostprocessBlock(PostprocessBlock &postprocessblock);

    void InitPPSOutput();
    void RunPostprocess(const double &time,const Mesh &mesh,const DofHandler &dofHandler,FE &fe,const SolutionSystem &solutionSystem);

    void CheckWhetherPPSIsValid(const Mesh &mesh);

    void PrintPostprocessInfo()const;


private:
    //******************************************************
    //*** for node pps
    //******************************************************
    double NodeValuePostProcess(const int &nodeid,string variablename,
                                const Mesh &mesh,const DofHandler &dofHandler,const SolutionSystem &solutionSystem);

    //******************************************************
    //*** for element pps
    //******************************************************
    double ElementValuePostProcess(const int &elmtid,string variablename,
                                   const Mesh &mesh,const DofHandler &dofHandler,const SolutionSystem &solutionSystem);
    //****************************************************************
    //*** do the elemental integration over specific domain name
    //****************************************************************
    double ElementalIntegralPostProcess(vector<string> domainnamelist,string variablename,
                                        const Mesh &mesh,const DofHandler &dofHandler,FE &fe,const SolutionSystem &solutionSystem);

    double AreaPostProcess(vector<string> sidenamelist,const Mesh &mesh,FE &fe);
    double VolumePostProcess(vector<string> domainnamelist,const Mesh &mesh,FE &fe);

    //****************************************************************
    //*** do the side integration over specific side-set name for
    //*** DoFs
    //****************************************************************
    double SideIntegralPostProcess(vector<string> sidenamelist,string dofname,
                                   const Mesh &mesh,const DofHandler &dofHandler,FE &fe,const SolutionSystem &solutionSystem);
    //****************************************************************
    //*** do the side integration over specific side-set name for
    //*** projected variable
    //****************************************************************
    double ProjVariableSideIntegralPostProcess(vector<string> sidenamelist,string variablename,
                                               const Mesh &mesh,FE &fe,const SolutionSystem &solutionSystem);

    //****************************************************************
    //*** do the side integration over specific side-set name for
    //*** projected variable
    //****************************************************************
    double Rank2MateSideIntegralPostProcess(vector<string> sidenamelist,string matename,const int &ii,const int &jj,
                                            const Mesh &mesh,FE &fe,const SolutionSystem &solutionSystem);

private:
    vector<PostprocessBlock> _PostProcessBlockList;
    int _nPostProcessBlocks;
    string _CSVFileName;
    string _InputFileName;
    vector<string> _VariableNameList;
    vector<double> _PPSValues;
    int _OutputInterval;

private:
    //*************************************************
    //*** for PETSc
    //*************************************************
    PetscMPIInt _rank;
    VecScatter _scatteru,_scatterproj;
    Vec _Useq,_ProjSeq;

};