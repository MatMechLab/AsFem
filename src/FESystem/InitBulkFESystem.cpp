//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.11.30
//+++ Purpose: Initialize all the FE system for bulk and interface
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/FESystem.h"

void FESystem::InitBulkFESystem(const Mesh &mesh,
                            const DofHandler &dofHandler,
                            FE &fe,
                            SolutionSystem &solution){

    _elNodes.InitNodes(mesh.GetBulkMeshNodesNumPerBulkElmt());

    _elConn.reserve(mesh.GetBulkMeshNodesNumPerBulkElmt());
    for(int i=0;i<mesh.GetBulkMeshNodesNumPerBulkElmt();++i){
        _elConn.push_back(0);
    }
    _elDofs.reserve(dofHandler.GetMaxDofsNumPerBulkElmt());
    _elDofsActiveFlag.reserve(dofHandler.GetMaxDofsNumPerBulkElmt());
    _elU.reserve(dofHandler.GetMaxDofsNumPerBulkElmt());
    _elV.reserve(dofHandler.GetMaxDofsNumPerBulkElmt());
    _elUold.reserve(dofHandler.GetMaxDofsNumPerBulkElmt());
    _elVold.reserve(dofHandler.GetMaxDofsNumPerBulkElmt());
    for(int i=0;i<dofHandler.GetMaxDofsNumPerBulkElmt();++i){
        _elDofs.push_back(0);
        _elDofsActiveFlag.push_back(1.0);
        _elU.push_back(0.0);
        _elV.push_back(0.0);
        _elUold.push_back(0.0);
        _elVold.push_back(0.0);
    }

    _gpU.reserve(dofHandler.GetDofsNumPerNode()+1);
    _gpV.reserve(dofHandler.GetDofsNumPerNode()+1);
    _gpUOld.reserve(dofHandler.GetDofsNumPerNode()+1);
    _gpVOld.reserve(dofHandler.GetDofsNumPerNode()+1);
    _gpGradU.reserve(dofHandler.GetDofsNumPerNode()+1);
    _gpGradV.reserve(dofHandler.GetDofsNumPerNode()+1);
    _gpGradUOld.reserve(dofHandler.GetDofsNumPerNode()+1);
    _gpGradVOld.reserve(dofHandler.GetDofsNumPerNode()+1);
    for(int i=0;i<dofHandler.GetDofsNumPerNode()+1;++i){
        _gpU.push_back(0.0);_gpUOld.push_back(0.0);
        _gpV.push_back(0.0);_gpVOld.push_back(0.0);
        _gpGradU.push_back(Vector3d(0.0));
        _gpGradV.push_back(Vector3d(0.0));
        _gpGradUOld.push_back(Vector3d(0.0));
        _gpGradVOld.push_back(Vector3d(0.0));
    }

    
    _nHist=solution.GetHistNumPerGPoint();
    _nProj=solution.GetProjNumPerNode();
    // the actual length of local gp's hist and proj can be much larger than the real one
    // but when we do the project, we only use the real one!!!
    _gpHist.reserve(20);
    _gpHistOld.reserve(20);
    for(int i=0;i<20;++i){
        _gpHist.push_back(0.0);
        _gpHistOld.push_back(0.0);
    }

    _gpProj.clear();
    

    _nGPoints=fe._BulkQPoint.GetQpPointsNum();
    
    
    _localK.Resize(dofHandler.GetMaxDofsNumPerBulkElmt(),dofHandler.GetMaxDofsNumPerBulkElmt());
    _localR.Resize(dofHandler.GetMaxDofsNumPerBulkElmt());

    _subK.Resize(dofHandler.GetDofsNumPerNode(),dofHandler.GetDofsNumPerNode());
    _subR.Resize(dofHandler.GetDofsNumPerNode());

    _K.resize(dofHandler.GetMaxDofsNumPerBulkElmt()*dofHandler.GetMaxDofsNumPerBulkElmt(),0.0);
    _R.resize(dofHandler.GetMaxDofsNumPerBulkElmt(),0.0);

    _localK.setZero();_localR.setZero();
    _subK.setZero();_subR.setZero();

    // set the factor to Ax=F system(this factor should be mesh dependent)
    // in order to get the most suitable one, we try to use 10 elements from the bulk
    int e,gpInd;
    double w,xi,eta,zeta,DetJac,JxW;
    int nDim=mesh.GetDim();
    _KMatrixFactor=1.0e16;
    int einc=int(1.0*mesh.GetBulkMeshBulkElmtsNum()/500);
    if(einc<1) einc=1;
    _BulkVolumes=0.0;
    for(e=1;e<=mesh.GetBulkMeshBulkElmtsNum();e+=einc){
        mesh.GetBulkMeshIthBulkElmtNodes(e,_elNodes);
        mesh.GetBulkMeshIthBulkElmtConn(e,_elConn);
        for(gpInd=1;gpInd<=fe._BulkQPoint.GetQpPointsNum();++gpInd){
            w=fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,0);
            xi=fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,1);
            if(nDim==1){
                fe._BulkShp.Calc(xi,_elNodes,true);
            }
            else if(nDim==2){
                eta=fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,2);
                fe._BulkShp.Calc(xi,eta,_elNodes,true);
            }
            else if(nDim==3){
                eta=fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,2);
                zeta=fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,3);
                fe._BulkShp.Calc(xi,eta,zeta,_elNodes,true);
            }
            DetJac=fe._BulkShp.GetDetJac();
            // JxW=1.0e3*DetJac*w; // it seems this is too small, it may lead the SNES solver failed
            //JxW=1.0e6*DetJac*w;
            JxW=DetJac*w;
            JxW=1.0; // maybe just keep it, do not scale the Ax=F system !!!
            if(JxW<_KMatrixFactor) _KMatrixFactor=JxW;
        }
    }
    
}