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
//+++ Date   : 2021.02.22
//+++ Purpose: run different types of postprocess for our FEM analysis
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocess.h"

void Postprocess::RunPostprocess(const double &time,const Mesh &mesh,const DofHandler &dofHandler,FE &fe,const SolutionSystem &solutionSystem){
    if(_nPostProcessBlocks<1){
        return;
    }
    PostprocessType ppstype;
    int nodeid,elmtid,iInd,jInd;
    string dofname,projvarname,rank2matename;
    vector<string> sidenamelist,domainnamelsit;
    for(int i=0;i<_nPostProcessBlocks;i++){
        ppstype=_PostProcessBlockList[i]._PostprocessType;
        nodeid=_PostProcessBlockList[i]._NodeID;
        elmtid=_PostProcessBlockList[i]._ElementID;
        iInd=_PostProcessBlockList[i]._iInd;
        jInd=_PostProcessBlockList[i]._jInd;
        sidenamelist=_PostProcessBlockList[i]._BoundaryNameList;
        domainnamelsit=_PostProcessBlockList[i]._DomainNameList;
        dofname=_PostProcessBlockList[i]._VariableName;
        projvarname=_PostProcessBlockList[i]._ProjVariableName;
        rank2matename=_PostProcessBlockList[i]._Rank2MateName;
        switch (ppstype) {
            case PostprocessType::NULLPPS:
                break;
            case PostprocessType::NODALVALUEPPS:
                _PPSValues[i]=NodeValuePostProcess(nodeid,dofname,mesh,dofHandler,solutionSystem);
                break;
            case PostprocessType::ELEMENTVALUEPPS:
                _PPSValues[i]=ElementValuePostProcess(elmtid,dofname,mesh,dofHandler,solutionSystem);
                break;
            case PostprocessType::AREAPPS:
                _PPSValues[i]=AreaPostProcess(sidenamelist,mesh,fe);
                break;
            case PostprocessType::SIDEINTEGRALPPS:
                _PPSValues[i]=SideIntegralPostProcess(sidenamelist,dofname,mesh,dofHandler,fe,solutionSystem);
                break;
            case PostprocessType::RANK2MATESIDEINTEGRALPPS:
                _PPSValues[i]=Rank2MateSideIntegralPostProcess(sidenamelist,rank2matename,iInd,jInd,mesh,fe,solutionSystem);
                break;
            case PostprocessType::ELEMENTINTEGRALPPS:
                _PPSValues[i]=ElementalIntegralPostProcess(domainnamelsit,dofname,mesh,dofHandler,fe,solutionSystem);
                break;
            case PostprocessType::VOLUMEPPS:
                _PPSValues[i]=VolumePostProcess(domainnamelsit,mesh,fe);
                break;
            case PostprocessType::PROJVARIABLESIDEINTEGRALPPS:
                _PPSValues[i]=ProjVariableSideIntegralPostProcess(sidenamelist,projvarname,mesh,fe,solutionSystem);
                break;
            default:
                MessagePrinter::PrintErrorTxt("unsupported postprocess type in RunPostprocess, please check your input file");
                MessagePrinter::AsFem_Exit();
                break;
        }
    }

    MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);
    if(_rank==0){
        ofstream out;
        out.open(_CSVFileName,ios::app|ios::out);
        if (!out.is_open()){
            string str="can\'t open the csv file(="+_CSVFileName+")!, please make sure you have write permission";
            MessagePrinter::PrintErrorTxt(str);
            MessagePrinter::AsFem_Exit();
        }
        char buff[20];
        sprintf(buff,"%-15.8e",time);
        out<<buff;
        for(const auto &it:_PPSValues){
            sprintf(buff,",%-15.8e",it);
            out<<buff;
        }
        out<<endl;
        out.close();
    }
}

