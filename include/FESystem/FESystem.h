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
//+++ Date   : 2020.11.29
//+++ Purpose: Implement the general tasks of FEM calculation in AsFem,
//+++          i.e. compute residual, compute jacobian
//+++          projection from gauss point to nodal point
//+++          assemble from local element to global, ...
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

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
#include "BCSystem/BCSystem.h"

#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"
#include "SolutionSystem/SolutionSystem.h"

#include "FE/FE.h"
#include "FE/ShapeFun.h"

#include "Utils/Vector3d.h"
#include "Utils/VectorXd.h"
#include "Utils/MatrixXd.h"

using namespace std;

class FESystem{
public:
    FESystem();
    void InitBulkFESystem(const Mesh &mesh,
                    const DofHandler &dofHandler,
                    FE &fe,
                    SolutionSystem &solution);

public:
    void SetHistNumPerGPoint(const int &n){_nHist=n;}
    void SetProjNumPerNode(const int &n) {_nProj=n;}
    void SetKMatrixFactor(const double &factor){_KMatrixFactor=factor;}
    inline double GetKMatrixFactor() const{return _KMatrixFactor;}
    void ResetMaxAMatrixValue(){_MaxKMatrixValue=-1.0e6;}
    void SetMaxAMatrixValue(const double &val) {_MaxKMatrixValue=val;}
    inline double GetMaxAMatrixValue()const {return _MaxKMatrixValue;}
    inline double GetBulkVolume() const {return _BulkVolumes;}

    // for FEM simulation related functions
    void FormBulkFE(const FECalcType &calctype,const double &t,const double &dt,const double (&ctan)[2],
                Mesh &mesh,const DofHandler &dofHandler,FE &fe,
                ElmtSystem &elmtSystem,MateSystem &mateSystem,
                const Vec &U,const Vec &V,
                Vec &Hist,Vec &HistOld,Vec &Proj,
                Mat &AMATRIX,Vec &RHS);
    
    
private:
    //*********************************************************
    //*** assemble residual to local and global one
    //*********************************************************
    void AssembleSubResidualToLocalResidual(const int &ndofspernode,const int &dofs,const int &iInd,
                                            const VectorXd &subR,VectorXd &localR);
    void AccumulateLocalResidual(const int &dofs,const vector<double> &dofsactiveflag,const double &JxW,
                                 const VectorXd &localR,vector<double> &sumR);
    void AssembleLocalResidualToGlobalResidual(const int &ndofs,const vector<int> &dofindex,
                                            const vector<double> &residual,Vec &rhs);

    //*********************************************************
    //*** assemble jacobian to local and global one
    //*********************************************************
    void AssembleSubJacobianToLocalJacobian(const int &ndofspernode,
                                            const int &iInd,const int &jInd,
                                            const MatrixXd &subK,MatrixXd &localK);
    void AccumulateLocalJacobian(const int &dofs,const vector<double> &dofsactiveflag,const double &JxW,
                                 const MatrixXd &localK,vector<double> &sumK);
    void AssembleLocalJacobianToGlobalJacobian(const int &ndofs,const vector<int> &dofindex,
                                            const vector<double> &jacobian,Mat &K);

    //*********************************************************
    //*** for projection
    //*********************************************************
    void AssembleLocalProjectionToGlobal(const int &nNodes,const double &DetJac,const ShapeFun &shp,const vector<double> &elProj,Vec &ProjVec);
    void Projection(const int &nTotalNodes,const int &nproj,Vec &ProjVec);


    void AssembleLocalHistToGlobal(const int &elmtid,const int &gpInd,const vector<double> &gpHist,Vec &Hist);
    
    void AssembleLocalToGlobal(const int &isw,const int &ndofs,vector<int> &elDofs,
                               vector<double> &localK,vector<double> &localR,
                               Mat &AMATRIX,Vec &RHS);

    //*********************************************************
    //*** for history variables
    //*********************************************************
    void AssembleLocalHistToGlobal(const int &e,const int &nhist,const int &ngp,const int &gpInd,const vector<double> &localHist,Vec &Hist);
    

public:
    void PrintFESystemInfo() const;



private:
    double _BulkVolumes=0.0;
    MatrixXd _localK;    //used in uel, the size is the maximum dofs per element
    VectorXd _localR;    //used in uel, the size is the maximum dofs per element
    MatrixXd _subK; // used in each sub element, the size is the maximum dofs per node
    VectorXd _subR; // used in each sub element, the size is the maximum dofs per node
    vector<double> _K,_R;//used in assemble
    
    
    Nodes _elNodes;
    vector<int> _elConn,_elDofs;
    vector<double> _elDofsActiveFlag;
    vector<double> _elU,_elV;
    vector<double> _gpU,_gpV;
    vector<double> _gpHist,_gpHistOld;
    map<string,double> _gpProj;
    vector<Vector3d> _gpGradU,_gpGradV;
    vector<double> _MaterialValues;
    Vector3d _gpCoord;
    int _nHist,_nProj,_nGPoints;
    double _MaxKMatrixValue=-1.0e9,_KMatrixFactor=0.1;

    ElmtType elmttype;
    MateType matetype;
    vector<int> localDofIndex;
    int mateindex;

private:
    //************************************
    //*** For PETSc related vairables
    PetscMPIInt _rank,_size;
    VecScatter _scatteru,_scatterv,_scatterproj,_scatterhist,_scatterhistold;
    Vec _Useq;// this can contain the ghost node from other processor
    Vec _Vseq;
    Vec _ProjSeq,_HistSeq,_HistOldSeq;
};