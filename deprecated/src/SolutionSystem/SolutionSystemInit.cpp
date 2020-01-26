#include "SolutionSystem/SolutionSystem.h"

void SolutionSystem::Init(){
    
    _UNorms.clear();
    _Utemp=Eigen::VectorXd(_nDofs);
    _Unew=Eigen::VectorXd(_nDofs);
    _U=Eigen::VectorXd(_nDofs);
    _Uold=Eigen::VectorXd(_nDofs);
    _Uolder=Eigen::VectorXd(_nDofs);
    _V=Eigen::VectorXd(_nDofs);
    _Vold=Eigen::VectorXd(_nDofs);
    _Volder=Eigen::VectorXd(_nDofs);
    
    _A=Eigen::VectorXd(_nDofs);
    _Aold=Eigen::VectorXd(_nDofs);
    _Aolder=Eigen::VectorXd(_nDofs);

    // must do this, otherwise Eigen will give your nan errors!!!
    _Utemp.setZero();_Unew.setZero();
    _U.setZero();_Uold.setZero();_Uolder.setZero();
    _V.setZero();_Vold.setZero();_Volder.setZero();
    

    if(_HasProjNameList){
        _ProjValue=Eigen::MatrixXd(_nNodes,1+_nProjsPerNode); // projected value of nodal point
    }
    else{
        _nProjsPerNode=9;
        _ProjValue=Eigen::MatrixXd(_nNodes,1+_nProjsPerNode);
        _ProjNameList.clear();
        _ProjNameList.reserve(_nProjsPerNode);
        string name;
        for(int i=1;i<=_nProjsPerNode;++i){
            name="Proj-"+to_string(i);
            _ProjNameList.push_back(name);
        }
    }
    _ProjValue.setZero();

    _nHistValuePerGPoint=5;
    _Hist=Eigen::MatrixXd(_nBulkElmts,_nHistValuePerGPoint*_nGPointPerBulkElmt);
    _HistOld=Eigen::MatrixXd(_nBulkElmts,_nHistValuePerGPoint*_nGPointPerBulkElmt);
    _Hist.setZero();
    _HistOld.setZero();
}