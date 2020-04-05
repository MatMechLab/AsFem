//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FESystem/FESystem.h"

void FESystem::InitFESystem(Mesh &mesh,
                            DofHandler &dofHandler,
                            FE &fe,
                            Solution &solution){

    _elNodes.InitNodes(mesh.GetNodesNumPerBulkElmt());

    _elConn.reserve(mesh.GetNodesNumPerBulkElmt());
    for(int i=0;i<mesh.GetNodesNumPerBulkElmt();++i){
        _elConn.push_back(0);
    }
    _elDofs.reserve(dofHandler.GetMaxDofsNumPerElmt());
    _elDofsActiveFlag.reserve(dofHandler.GetMaxDofsNumPerElmt());
    _elU.reserve(dofHandler.GetMaxDofsNumPerElmt());
    _elV.reserve(dofHandler.GetMaxDofsNumPerElmt());
    for(int i=0;i<dofHandler.GetMaxDofsNumPerElmt();++i){
        _elDofs.push_back(0);
        _elDofsActiveFlag.push_back(1.0);
        _elU.push_back(0.0);
        _elV.push_back(0.0);
    }

    _gpU.reserve(dofHandler.GetDofsNumPerNode());
    _gpV.reserve(dofHandler.GetDofsNumPerNode());
    _gpGradU.reserve(dofHandler.GetDofsNumPerNode());
    _gpGradV.reserve(dofHandler.GetDofsNumPerNode());
    for(int i=0;i<dofHandler.GetDofsNumPerNode();++i){
        _gpU.push_back(0.0);
        _gpV.push_back(0.0);
        _gpGradU.push_back(Vector3d(0.0));
        _gpGradV.push_back(Vector3d(0.0));
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
    _gpProj.reserve(20);
    for(int i=0;i<20;++i){
        _gpProj.push_back(0.0);
    }

    _nGPoints=fe._qp_bulk.GetQpPointsNum();
    
    
    

    
    _localK=MatrixXd(dofHandler.GetMaxDofsNumPerElmt(),dofHandler.GetMaxDofsNumPerElmt());
    _localR=VectorXd(dofHandler.GetMaxDofsNumPerElmt());

    _K.resize(dofHandler.GetMaxDofsNumPerElmt()*dofHandler.GetMaxDofsNumPerElmt(),0.0);
    _R.resize(dofHandler.GetMaxDofsNumPerElmt(),0.0);

    _localK.setZero();_localR.setZero();

    // set the factor to Ax=F system(this factor should be mesh dependent)
    // in order to get the most suitable one, we try to use 10 elements from the bulk
    int e,gpInd;
    double w,xi,eta,zeta,DetJac,JxW;
    int nDim=mesh.GetDim();
    _KMatrixFactor=1.0e16;
    int einc=int(1.0*mesh.GetBCElmtsNum()/10);
    if(einc<1) einc=1;
    for(e=1;e<=mesh.GetBulkElmtsNum();e+=einc){
        mesh.GetIthBulkElmtNodes(e,_elNodes);
        mesh.GetIthBulkElmtConn(e,_elConn);
        for(gpInd=1;gpInd<=fe._qp_bulk.GetQpPointsNum();++gpInd){
            w=fe._qp_bulk.GetIthQpPointJthCoord(gpInd,0);
            xi=fe._qp_bulk.GetIthQpPointJthCoord(gpInd,1);
            if(nDim==1){
                fe._shp_bulk.Calc(xi,_elNodes,true);
            }
            else if(nDim==2){
                eta=fe._qp_bulk.GetIthQpPointJthCoord(gpInd,2);
                fe._shp_bulk.Calc(xi,eta,_elNodes,true);
            }
            else if(nDim==3){
                eta=fe._qp_bulk.GetIthQpPointJthCoord(gpInd,2);
                zeta=fe._qp_bulk.GetIthQpPointJthCoord(gpInd,3);
                fe._shp_bulk.Calc(xi,eta,zeta,_elNodes,true);
            }
            DetJac=fe._shp_bulk.GetDetJac();
            // JxW=1.0e3*DetJac*w; // it seems this is too small, it may lead the SNES solver failed
            //JxW=1.0e6*DetJac*w;
            JxW=DetJac*w;
            JxW=1.0; // maybe just keep it, do not scale the Ax=F system !!!
            if(JxW<_KMatrixFactor) _KMatrixFactor=JxW;
        }
    }
    
    
}