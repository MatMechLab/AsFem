//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_FESYSTEM_H
#define ASFEM_FESYSTEM_H


#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>

#include "petsc.h"


//*************************************
//*** For AsFem's own header file
//*************************************
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "BCs/BCSystem.h"

#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"
#include "Solution/Solution.h"

#include "FE/FE.h"

#include "Utils/Vector3d.h"
#include "Utils/VectorXd.h"
#include "Utils/MatrixXd.h"

using namespace std;

class FESystem{
public:
    FESystem();
    void InitFESystem(Mesh &mesh,
                      DofHandler &dofHandler,
                      FE &fe,
                      Solution &solution);

public:
    void SetHistNumPerGPoint(const int &n){_nHist=n;}
    void SetProjNumPerNode(const int &n) {_nProj=n;}
    void SetKMatrixFactor(const double &factor){_KMatrixFactor=factor;}
    inline double GetKMatrixFactor() const{return _KMatrixFactor;}
    void ResetMaxAMatrixValue(){_MaxKMatrixValue=-1.0e6;}
    void SetMaxAMatrixValue(const double &val) {_MaxKMatrixValue=val;}
    inline double GetMaxAMatrixValue()const {return _MaxKMatrixValue;}
    inline double GetVolume() const {return _Volumes;}
    // for FEM simulation related functions
    void FormFE(const int &isw,const double &t,const double &dt,const double (&ctan)[2],
                Mesh &mesh,DofHandler &dofHandler,FE &fe,
                ElmtSystem &elmtSystem,MateSystem &mateSystem,
                const Vec &U,const Vec &V,
                Vec &Hist,const Vec &HistOld,Vec &Proj,
                Mat &AMATRIX,Vec &RHS);
    
    
private:
    void Projection(const int &nNodes,Vec &Proj,const vector<double> &elProj,ShapeFun &shp,const double &DetJac);
    void AssembleLocalHistToGlobal(const int &elmtid,const int &gpInd,const vector<double> &gpHist,Vec &Hist);
    
    void AssembleLocalToGlobal(const int &isw,const int &ndofs,vector<int> &elDofs,
                               vector<double> &localK,vector<double> &localR,
                               Mat &AMATRIX,Vec &RHS);
    

public:
    void PrintFESystemInfo() const;



private:
    double _Volumes=0.0;
    MatrixXd _localK;    //used in uel
    VectorXd _localR;    //used in uel
    vector<double> _K,_R;//used in assemble
    
    
    Nodes _elNodes;
    vector<int> _elConn,_elDofs;
    vector<double> _elDofsActiveFlag;
    vector<double> _elU,_elV;
    vector<double> _gpU,_gpV;
    vector<double> _gpHist,_gpHistOld,_gpProj;
    vector<Vector3d> _gpGradU,_gpGradV;
    vector<double> _MaterialValues;
    Vector3d _gpCoord;
    int _nHist,_nProj,_nGPoints;
    double _MaxKMatrixValue=-1.0e3,_KMatrixFactor=0.1;

private:
    //************************************
    //*** For PETSc related vairables
    PetscMPIInt _rank,_size;
    VecScatter _scatteru,_scatterv,_scatterproj,_scatterhist,_scatterhistold;
    Vec _Useq;// this can contain the ghost node from other processor
    Vec _Vseq;
    Vec _ProjSeq,_HistSeq,_HistOldSeq;
};


#endif // ASFEM_FESYSTEM_H