//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.07.22
//+++ Purpose: Implement least square projection from integration
//+++          points to nodal points
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ProjectionSystem/LeastSquareProjection.h"

LeastSquareProjection::LeastSquareProjection(){

}
void LeastSquareProjection::initMyProjection(const FECell &t_FECell, const DofHandler &t_DofHandler) {
    /**
     * for the qpoints used in projection system
     */
    int m_max_nodal_dofs=t_DofHandler.getMaxDofsPerNode();
    m_BulkElmtNodesNum=t_FECell.getFECellNodesNumPerBulkElmt();

    m_ElmtConn.resize(m_BulkElmtNodesNum,0);
    m_SubElmtDofIDs.resize(m_max_nodal_dofs+1,0);

    m_LocalElmtSoln.m_QpU.resize(m_max_nodal_dofs+1,0.0);
    m_LocalElmtSoln.m_QpUold.resize(m_max_nodal_dofs+1,0.0);
    m_LocalElmtSoln.m_QpUolder.resize(m_max_nodal_dofs+1,0.0);

    m_LocalElmtSoln.m_QpV.resize(m_max_nodal_dofs+1,0.0);
    m_LocalElmtSoln.m_QpA.resize(m_max_nodal_dofs+1,0.0);

    m_LocalElmtSoln.m_QpGradU.resize(m_max_nodal_dofs+1,Vector3d(0.0));
    m_LocalElmtSoln.m_QpGradUold.resize(m_max_nodal_dofs+1,Vector3d(0.0));
    m_LocalElmtSoln.m_QpGradUolder.resize(m_max_nodal_dofs+1,Vector3d(0.0));

    m_LocalElmtSoln.m_QpGradV.resize(m_max_nodal_dofs+1,Vector3d(0.0));

    m_LocalElmtSoln.m_Qpgradu.resize(m_max_nodal_dofs+1,Vector3d(0.0));
    m_LocalElmtSoln.m_Qpgradv.resize(m_max_nodal_dofs+1,Vector3d(0.0));

    m_LocalElmtInfo.m_QpCoords=0.0;
    m_LocalElmtInfo.m_QpCoords0=0.0;

    m_Nodes.resize(m_BulkElmtNodesNum);
    m_Nodes0.resize(m_BulkElmtNodesNum);
}

//******************************************************
//*** for local projection action
//******************************************************
void LeastSquareProjection::localProjectionAction(const int &NodesNum,
                                                  const vector<int> &ElConn,
                                                  const double &DetJac,
                                                  const ShapeFun &Shp,
                                                  const MaterialsContainer &Mate,
                                                  ProjectionData &Data){
    
    if(Data.m_ScalarProjMateNum){
        projectLocalScalarMate2Global(NodesNum,ElConn,DetJac,Shp,Mate,Data.m_ScalarProjMateNum,Data.m_ScalarProjMateNameList,Data.m_ScalarProjMateVecList);
    }
    if(Data.m_VectorProjMateNum){
        projectLocalVectorMate2Global(NodesNum,ElConn,DetJac,Shp,Mate,Data.m_VectorProjMateNum,Data.m_VectorProjMateNamelist,Data.m_VectorProjMateVecList);
    }
    if(Data.m_Rank2ProjMateNum){
        projectLocalRank2Mate2Global(NodesNum,ElConn,DetJac,Shp,Mate,Data.m_Rank2ProjMateNum,Data.m_Rank2ProjMateNameList,Data.m_Rank2ProjMateVecList);
    }
    if(Data.m_Rank4ProjMateNum){
        projectLocalRank4Mate2Global(NodesNum,ElConn,DetJac,Shp,Mate,Data.m_Rank4ProjMateNum,Data.m_Rank4ProjMateNameList,Data.m_Rank4ProjMateVecList);
    }
}
//******************************************************
void LeastSquareProjection::projectLocalScalarMate2Global(const int &NodesNum,
                                                          const vector<int> &ElConn,
                                                          const double &DetJac,
                                                          const ShapeFun &Shp,
                                                          const MaterialsContainer &Mate,
                                                          const int &ProjNum,
                                                          const vector<string> &ScalarMateNameList,
                                                          vector<Vector> &ScalarProjMateVecList){
    if(ProjNum!=static_cast<int>(ScalarMateNameList.size())){
        MessagePrinter::printErrorTxt("the number of projected scalar material does not match with the given one, please check your code");
        MessagePrinter::exitAsFem();
    }
    bool HasName;
    double w,mateval;
    int j,k,jInd,iInd;
    const auto &scalmat=Mate.getScalarMaterialsRef();// const ref: no per-node map copy
    for(k=1;k<=ProjNum;k++){
        // locate the projected material value ONCE (was copied+searched for every node)
        HasName=false; mateval=0.0;
        for(const auto &it:scalmat){
            if(it.first==ScalarMateNameList[k-1]){ HasName=true; mateval=it.second; break; }
        }
        if(!HasName){
            MessagePrinter::printErrorTxt("can\'t find the scalar material ("+ScalarMateNameList[k-1]+") for projection, please check either your input file or your material code");
            MessagePrinter::exitAsFem();
        }
        for(j=1;j<=NodesNum;j++){
            iInd=ElConn[j-1];
            w=DetJac*Shp.shape_value(j);
            jInd=(iInd-1)*(1+1)+1;
            ScalarProjMateVecList[k-1].addValue(jInd,w);
            jInd=(iInd-1)*(1+1)+2;
            ScalarProjMateVecList[k-1].addValue(jInd,w*mateval);
        }
    }
}
//******************************************************
void LeastSquareProjection::projectLocalVectorMate2Global(const int &NodesNum,
                                                          const vector<int> &ElConn,
                                                          const double &DetJac,
                                                          const ShapeFun &Shp,
                                                          const MaterialsContainer &Mate,
                                                          const int &ProjNum,
                                                          const vector<string> &VectorMateNameList,
                                                          vector<Vector> &VectorProjMateVecList){
    if(ProjNum!=static_cast<int>(VectorMateNameList.size())){
        MessagePrinter::printErrorTxt("the number of projected vector material does not match with the given one, please check your code");
        MessagePrinter::exitAsFem();
    }
    bool HasName;
    double w;
    int j,k,jInd,iInd;
    Vector3d matevec;
    const auto &vecmat=Mate.getVectorMaterialsRef();// const ref: no per-node map copy
    for(k=1;k<=ProjNum;k++){
        HasName=false; matevec=0.0;
        for(const auto &it:vecmat){
            if(it.first==VectorMateNameList[k-1]){ HasName=true; matevec=it.second; break; }
        }
        if(!HasName){
            MessagePrinter::printErrorTxt("can\'t find the vector material ("+VectorMateNameList[k-1]+") for projection, please check either your input file or your material code");
            MessagePrinter::exitAsFem();
        }
        for(j=1;j<=NodesNum;j++){
            iInd=ElConn[j-1];
            w=DetJac*Shp.shape_value(j);
            jInd=(iInd-1)*(1+3)+1;
            VectorProjMateVecList[k-1].addValue(jInd,w);
            jInd=(iInd-1)*(1+3)+2;
            VectorProjMateVecList[k-1].addValue(jInd,matevec(1)*w);
            jInd=(iInd-1)*(1+3)+3;
            VectorProjMateVecList[k-1].addValue(jInd,matevec(2)*w);
            jInd=(iInd-1)*(1+3)+4;
            VectorProjMateVecList[k-1].addValue(jInd,matevec(3)*w);
        }
    }
}
//******************************************************
void LeastSquareProjection::projectLocalRank2Mate2Global(const int &NodesNum,
                                                         const vector<int> &ElConn,
                                                         const double &DetJac,
                                                         const ShapeFun &Shp,
                                                         const MaterialsContainer &Mate,
                                                         const int &ProjNum,
                                                         const vector<string> &Rank2MateNameList,
                                                         vector<Vector> &Rank2ProjMateVecList){
    if(ProjNum!=static_cast<int>(Rank2MateNameList.size())){
        MessagePrinter::printErrorTxt("the number of projected rank-2 material does not match with the given one, please check your code");
        MessagePrinter::exitAsFem();
    }
    bool HasName;
    double w;
    int j,k,jInd,iInd;
    int i1,j1,ii;
    Rank2Tensor mate2;
    const auto &rank2mat=Mate.getRank2MaterialsRef();// const ref: no per-node map copy
    for(k=1;k<=ProjNum;k++){
        HasName=false; mate2.setToZeros();
        for(const auto &it:rank2mat){
            if(it.first==Rank2MateNameList[k-1]){ HasName=true; mate2=it.second; break; }
        }
        if(!HasName){
            MessagePrinter::printErrorTxt("can\'t find the rank-2 material ("+Rank2MateNameList[k-1]+") for projection, please check either your input file or your material code");
            MessagePrinter::exitAsFem();
        }
        for(j=1;j<=NodesNum;j++){
            iInd=ElConn[j-1];
            w=DetJac*Shp.shape_value(j);
            jInd=(iInd-1)*(9+1)+0+1;
            Rank2ProjMateVecList[k-1].addValue(jInd,w);
            ii=0;
            for(i1=1;i1<=3;i1++){
                for(j1=1;j1<=3;j1++){
                    ii+=1;
                    jInd=(iInd-1)*(9+1)+ii+1;
                    Rank2ProjMateVecList[k-1].addValue(jInd,mate2(i1,j1)*w);
                }
            }
        }
    }
}
//******************************************************
void LeastSquareProjection::projectLocalRank4Mate2Global(const int &NodesNum,
                                                         const vector<int> &ElConn,
                                                         const double &DetJac,
                                                         const ShapeFun &Shp,
                                                         const MaterialsContainer &Mate,
                                                         const int &ProjNum,
                                                         const vector<string> &Rank4MateNameList,
                                                         vector<Vector> &Rank4ProjMateVecList){
    if(ProjNum!=static_cast<int>(Rank4MateNameList.size())){
        MessagePrinter::printErrorTxt("the number of projected rank-4 material does not match with the given one, please check your code");
        MessagePrinter::exitAsFem();
    }
    bool HasName;
    double w;
    int j,k,jInd,iInd;
    int i1,j1,ii;
    Rank4Tensor mate4;
    const auto &rank4mat=Mate.getRank4MaterialsRef();// const ref: no per-node map copy
    for(k=1;k<=ProjNum;k++){
        HasName=false; mate4.setToZeros();
        for(const auto &it:rank4mat){
            if(it.first==Rank4MateNameList[k-1]){ HasName=true; mate4=it.second; break; }
        }
        if(!HasName){
            MessagePrinter::printErrorTxt("can\'t find the rank-4 material ("+Rank4MateNameList[k-1]+") for projection, please check either your input file or your material code");
            MessagePrinter::exitAsFem();
        }
        for(j=1;j<=NodesNum;j++){
            iInd=ElConn[j-1];
            w=DetJac*Shp.shape_value(j);
            ii=0;
            jInd=(iInd-1)*(36+1)+ii+1;
            Rank4ProjMateVecList[k-1].addValue(jInd,w);
            for(i1=1;i1<=6;i1++){
                for(j1=1;j1<=6;j1++){
                    ii+=1;
                    jInd=(iInd-1)*(36+1)+ii+1;
                    Rank4ProjMateVecList[k-1].addValue(jInd,mate4.getVoigtComponent(i1,j1)*w);
                }
            }
        }
    }
}
//******************************************************
//*** for global projection action
//******************************************************
void LeastSquareProjection::globalProjectionAction(const FECell &t_FECell,
                                                   ProjectionData &Data){
    int j,k,iInd,nproj;
    double value,weight,newvalue;

    // assemble all the local projection value
    for (auto &it:Data.m_ScalarProjMateVecList) it.assemble();
    for (auto &it:Data.m_VectorProjMateVecList) it.assemble();
    for (auto &it:Data.m_Rank2ProjMateVecList) it.assemble();
    for (auto &it:Data.m_Rank4ProjMateVecList) it.assemble();

    // make a ghost copy for all of them
    for (auto &it:Data.m_ScalarProjMateVecList) it.makeGhostCopy();
    for (auto &it:Data.m_VectorProjMateVecList) it.makeGhostCopy();
    for (auto &it:Data.m_Rank2ProjMateVecList) it.makeGhostCopy();
    for (auto &it:Data.m_Rank4ProjMateVecList) it.makeGhostCopy();

    // set all of them to zero vectors
    for (auto &it:Data.m_ScalarProjMateVecList) it.setToZero();
    for (auto &it:Data.m_VectorProjMateVecList) it.setToZero();
    for (auto &it:Data.m_Rank2ProjMateVecList) it.setToZero();
    for (auto &it:Data.m_Rank4ProjMateVecList) it.setToZero();

    // read-only ghost array pointer for every projection vector: arr[iInd-1] replaces a per-value
    // VecGetValues for the weight/value reads below (the writes still go through addValue).
    vector<const PetscScalar*> _gS(Data.m_ScalarProjMateNum),_gV(Data.m_VectorProjMateNum),
                               _gR2(Data.m_Rank2ProjMateNum),_gR4(Data.m_Rank4ProjMateNum);
    for(int m=0;m<Data.m_ScalarProjMateNum;m++) _gS[m]=Data.m_ScalarProjMateVecList[m].getGhostArrayRead();
    for(int m=0;m<Data.m_VectorProjMateNum;m++) _gV[m]=Data.m_VectorProjMateVecList[m].getGhostArrayRead();
    for(int m=0;m<Data.m_Rank2ProjMateNum;m++) _gR2[m]=Data.m_Rank2ProjMateVecList[m].getGhostArrayRead();
    for(int m=0;m<Data.m_Rank4ProjMateNum;m++) _gR4[m]=Data.m_Rank4ProjMateVecList[m].getGhostArrayRead();

    for (const auto &GlobalNodeID:t_FECell.getLocalFECellNodeIDsCopy()) {
        //****************************************
        // for scalar materials
        //****************************************
        nproj = Data.m_ScalarProjMateNum;
        if (nproj) {
            for (j = 1; j <= nproj; j++) {
                iInd = (GlobalNodeID - 1) * (1 + 1) + 1;
                weight = _gS[j-1][iInd-1];
                iInd = (GlobalNodeID - 1) * (1 + 1) + 2;
                value = _gS[j-1][iInd-1];
                if (abs(weight) > 1.0e-15) {
                    newvalue = value / weight;
                } else {
                    newvalue = value /1.0e-8;
                }
                Data.m_ScalarProjMateVecList[j-1].addValue(iInd, newvalue);
            }
        }
        //****************************************
        // for vector materials
        //****************************************
        nproj = Data.m_VectorProjMateNum;
        if (nproj) {
            for (j = 1; j <= nproj; j++) {
                iInd = (GlobalNodeID - 1) * (1 + 3) + 0 + 1;
                weight = _gV[j-1][iInd-1];
                for (k = 1; k <= 3; k++) {
                    iInd = (GlobalNodeID - 1) * (3 + 1)+ k + 1;
                    value = _gV[j-1][iInd-1];
                    if (abs(weight) > 1.0e-15) {
                        newvalue = value / weight;
                    } else {
                        newvalue = value /1.0e-8;
                    }
                    Data.m_VectorProjMateVecList[j-1].addValue(iInd, newvalue);
                }
            }
        }
        //****************************************
        // for rank-2 materials
        //****************************************
        nproj = Data.m_Rank2ProjMateNum;
        if (nproj) {
            for (j = 1; j <= nproj; j++) {
                iInd = (GlobalNodeID - 1) * (1 + 9) + 0 + 1;
                weight = _gR2[j-1][iInd-1];
                for (k = 1; k <= 9; k++) {
                    iInd = (GlobalNodeID - 1) * (9 + 1) + k + 1;
                    value = _gR2[j-1][iInd-1];
                    if (abs(weight) > 1.0e-15) {
                        newvalue = value / weight;
                    } else {
                        newvalue = value /1.0e-8;
                    }
                    Data.m_Rank2ProjMateVecList[j-1].addValue(iInd, newvalue);
                }
            }
        }
        //****************************************
        // for rank-4 materials
        //****************************************
        nproj = Data.m_Rank4ProjMateNum;
        if (nproj) {
            for (j = 1; j <= nproj; j++) {
                iInd = (GlobalNodeID - 1) * (1 + 36) + 0 + 1;
                weight = _gR4[j-1][iInd-1];
                for (k = 1; k <= 36; k++) {
                    iInd = (GlobalNodeID - 1) * (36 + 1)+ k + 1;
                    value = _gR4[j-1][iInd-1];
                    if (abs(weight) > 1.0e-15) {
                        newvalue = value / weight;
                    } else {
                        newvalue = value /1.0e-8;
                    }
                    Data.m_Rank4ProjMateVecList[j-1].addValue(iInd, newvalue);
                }
            }
        }
    } //end-of-local-element-nodes-loop

    // release the ghost array pointers
    for(int m=0;m<Data.m_ScalarProjMateNum;m++) Data.m_ScalarProjMateVecList[m].restoreGhostArrayRead(&_gS[m]);
    for(int m=0;m<Data.m_VectorProjMateNum;m++) Data.m_VectorProjMateVecList[m].restoreGhostArrayRead(&_gV[m]);
    for(int m=0;m<Data.m_Rank2ProjMateNum;m++) Data.m_Rank2ProjMateVecList[m].restoreGhostArrayRead(&_gR2[m]);
    for(int m=0;m<Data.m_Rank4ProjMateNum;m++) Data.m_Rank4ProjMateVecList[m].restoreGhostArrayRead(&_gR4[m]);

    // finish the final assemble
    for (auto &it:Data.m_ScalarProjMateVecList) it.assemble();
    for (auto &it:Data.m_VectorProjMateVecList) it.assemble();
    for (auto &it:Data.m_Rank2ProjMateVecList) it.assemble();
    for (auto &it:Data.m_Rank4ProjMateVecList) it.assemble();

    // destroy the ghost copy
    for (auto &it:Data.m_ScalarProjMateVecList) it.destroyGhostCopy();
    for (auto &it:Data.m_VectorProjMateVecList) it.destroyGhostCopy();
    for (auto &it:Data.m_Rank2ProjMateVecList) it.destroyGhostCopy();
    for (auto &it:Data.m_Rank4ProjMateVecList) it.destroyGhostCopy();

}
void LeastSquareProjection::executeMyProjection(const FECell &t_FECell,
                                                const DofHandler &t_DofHandler,
                                                const ElmtSystem &t_ElmtSystem,
                                                MateSystem &t_MateSystem,
                                                FE &t_FE,
                                                SolutionSystem &t_SolnSystem,
                                                const FEControlInfo &t_FECtrlInfo,
                                                ProjectionData &Data) {
    // if no projection is required, then return back
    if(Data.m_ScalarProjMateNum<1 &&
       Data.m_VectorProjMateNum<1 &&
       Data.m_Rank2ProjMateNum<1 &&
       Data.m_Rank4ProjMateNum<1){
        return;
    }

    // before we do the projection, set all the vector to zeros
    for (auto &it:Data.m_ScalarProjMateVecList) it.setToZero();
    for (auto &it:Data.m_VectorProjMateVecList) it.setToZero();
    for (auto &it:Data.m_Rank2ProjMateVecList) it.setToZero();
    for (auto &it:Data.m_Rank4ProjMateVecList) it.setToZero();

    // for the current and the previous steps' solution array
    t_SolnSystem.m_Ucurrent.makeGhostCopy();//
    t_SolnSystem.m_Uold.makeGhostCopy();
    t_SolnSystem.m_Uolder.makeGhostCopy();

    // for the velocity and acceleration
    t_SolnSystem.m_V.makeGhostCopy();
    t_SolnSystem.m_A.makeGhostCopy();

    int QpointsNum;
    double xi,eta,zeta,w,J,JxW;

    int SubElmtBlockID;

    m_LocalElmtInfo.m_Dim=t_FECell.getFECellMaxDim();
    m_LocalElmtInfo.m_NodesNum=t_FECell.getFECellNodesNumPerBulkElmt();

    // reference (not a copy) of the local cells: the loop only reads them
    const vector<SingleMeshCell> &MyLocalCellVec=t_FECell.getLocalBulkFECellVecRef();
    m_BulkElmtNodesNum=t_FECell.getFECellNodesNumPerBulkElmt();

    QpointsNum=t_FE.m_BulkQpoints.getQPointsNum();
    m_LocalElmtInfo.m_QpointsNum=QpointsNum;

    // per-element cache of the nodal solution values, indexed by (node, dof); fetched ONCE per
    // element instead of from the ghost vectors at every quadrature point (a VecGetValues per access)
    const int _maxDofsPerNode=t_DofHandler.getMaxDofsPerNode();
    const int _cacheW=m_BulkElmtNodesNum*_maxDofsPerNode;
    vector<double> _cUolder(_cacheW,0.0),_cUold(_cacheW,0.0),_cU(_cacheW,0.0),_cV(_cacheW,0.0),_cA(_cacheW,0.0);

    // read-only pointers to the to-all ghost arrays; arr[gdof-1] replaces a per-value VecGetValues
    const PetscScalar *_aUolder=t_SolnSystem.m_Uolder.getGhostArrayRead();
    const PetscScalar *_aUold  =t_SolnSystem.m_Uold.getGhostArrayRead();
    const PetscScalar *_aU     =t_SolnSystem.m_Ucurrent.getGhostArrayRead();
    const PetscScalar *_aV     =t_SolnSystem.m_V.getGhostArrayRead();
    const PetscScalar *_aA     =t_SolnSystem.m_A.getGhostArrayRead();

    for (int e=1;e<=t_FECell.getLocalFECellBulkElmtsNum();e++) {
        m_LocalElmtInfo.m_Dim=MyLocalCellVec[e-1].Dim;
        m_LocalElmtInfo.m_NodesNum=MyLocalCellVec[e-1].NodesNumPerElmt;
        m_LocalElmtInfo.m_Dt=t_FECtrlInfo.Dt;
        m_LocalElmtInfo.m_T=t_FECtrlInfo.T;

        m_Nodes=MyLocalCellVec[e-1].ElmtNodeCoords;
        m_Nodes0=MyLocalCellVec[e-1].ElmtNodeCoords0;
        m_ElmtConn=MyLocalCellVec[e-1].ElmtConn;

        // build the per-element nodal-value cache once
        for(int jn=1;jn<=m_BulkElmtNodesNum;jn++){
            const int gnode=MyLocalCellVec[e-1].ElmtConn[jn-1];
            for(int d=1;d<=_maxDofsPerNode;d++){
                const int gdof=t_DofHandler.getIthNodeJthDofID(gnode,d);
                if(gdof>=1){
                    const int idx=(jn-1)*_maxDofsPerNode+(d-1);
                    _cUolder[idx]=_aUolder[gdof-1];
                    _cUold[idx]  =_aUold[gdof-1];
                    _cU[idx]     =_aU[gdof-1];
                    _cV[idx]     =_aV[gdof-1];
                    _cA[idx]     =_aA[gdof-1];
                }
            }
        }

        for (int qp=1;qp<=QpointsNum;qp++) {
            w=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,0);
            xi=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,1);
            if (m_LocalElmtInfo.m_Dim==1) {
                eta=0.0;zeta=0.0;
            }
            else if (m_LocalElmtInfo.m_Dim==2) {
                eta=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,2);
                zeta=0.0;
            }
            else if (m_LocalElmtInfo.m_Dim==3) {
                eta=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,2);
                zeta=t_FE.m_BulkQpoints.getIthPointJthCoord(qp,3);
            }
            t_FE.m_BulkShp.calc(xi,eta,zeta,m_Nodes0,true);
            J=t_FE.m_BulkShp.getJacDet();
            JxW=J*w;

            m_LocalElmtInfo.m_QpointID=qp;
            m_LocalElmtInfo.m_QpCoords0=0.0;
            for (int i=1;i<=m_LocalElmtInfo.m_NodesNum;i++) {
                m_LocalElmtInfo.m_QpCoords0(1)+=t_FE.m_BulkShp.shape_value(i)*m_Nodes0(i,1);
                m_LocalElmtInfo.m_QpCoords0(2)+=t_FE.m_BulkShp.shape_value(i)*m_Nodes0(i,2);
                m_LocalElmtInfo.m_QpCoords0(3)+=t_FE.m_BulkShp.shape_value(i)*m_Nodes0(i,3);
            }

            t_MateSystem.m_MaterialContainerOld.getScalarMaterialsRef()=t_SolnSystem.m_QpointsScalarMaterialsOld_Local[(e-1)*QpointsNum+qp-1];
            t_MateSystem.m_MaterialContainerOld.getVectorMaterialsRef()=t_SolnSystem.m_QpointsVectorMaterialsOld_Local[(e-1)*QpointsNum+qp-1];
            t_MateSystem.m_MaterialContainerOld.getRank2MaterialsRef()=t_SolnSystem.m_QpointsRank2MaterialsOld_Local[(e-1)*QpointsNum+qp-1];
            t_MateSystem.m_MaterialContainerOld.getRank4MaterialsRef()=t_SolnSystem.m_QpointsRank4MaterialsOld_Local[(e-1)*QpointsNum+qp-1];

            for (int SubElmt=1;SubElmt<=t_ElmtSystem.getLocalIthBulkElmtSubElmtsNum(e);SubElmt++) {
                SubElmtBlockID=t_ElmtSystem.getLocalIthBulkElmtJthSubElmtID(e,SubElmt);
                m_SubElmtDofs=static_cast<int>(t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_DofIDs.size());
                m_LocalElmtInfo.m_DofsNum=m_SubElmtDofs;

                //*********************************************************************
                //*** calculate physical quantities on each integration point
                //*********************************************************************
                // here, the U/V/A vector's index should start from 1, to make it consistent with
                // 1st dof, 2nd dof, 3rd one , and so on. So, please let the index starts from 1 !!!
                for (int i=1;i<=m_LocalElmtInfo.m_DofsNum;i++){
                    m_SubElmtDofIDs[i-1]=t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_DofIDs[i-1];

                    m_LocalElmtSoln.m_QpU[i]=0.0;
                    m_LocalElmtSoln.m_QpUold[i]=0.0;
                    m_LocalElmtSoln.m_QpUolder[i]=0.0;

                    m_LocalElmtSoln.m_QpV[i]=0.0;
                    m_LocalElmtSoln.m_QpA[i]=0.0;

                    m_LocalElmtSoln.m_QpGradU[i]=0.0;
                    m_LocalElmtSoln.m_QpGradUold[i]=0.0;
                    m_LocalElmtSoln.m_QpGradUolder[i]=0.0;

                    m_LocalElmtSoln.m_QpGradV[i]=0.0;
                    for (int j=1;j<=m_BulkElmtNodesNum;j++) {
                        // read the (node j, sub-dof i) nodal values from the per-element cache
                        // (each value fetched once per element instead of re-read every qpoint,
                        //  and once per value instead of once per gradient component)
                        const int cidx=(j-1)*_maxDofsPerNode+(m_SubElmtDofIDs[i-1]-1);
                        const double uolder=_cUolder[cidx];
                        const double uold  =_cUold[cidx];
                        const double u     =_cU[cidx];
                        const double v     =_cV[cidx];
                        const double a     =_cA[cidx];
                        const double N=t_FE.m_BulkShp.shape_value(j);
                        const Vector3d dN=t_FE.m_BulkShp.shape_grad(j);

                        m_LocalElmtSoln.m_QpUolder[i]+=N*uolder;
                        m_LocalElmtSoln.m_QpUold[i]  +=N*uold;
                        m_LocalElmtSoln.m_QpU[i]     +=N*u;
                        m_LocalElmtSoln.m_QpV[i]     +=N*v;
                        m_LocalElmtSoln.m_QpA[i]     +=N*a;

                        m_LocalElmtSoln.m_QpGradUolder[i](1)+=dN(1)*uolder;
                        m_LocalElmtSoln.m_QpGradUolder[i](2)+=dN(2)*uolder;
                        m_LocalElmtSoln.m_QpGradUolder[i](3)+=dN(3)*uolder;
                        //
                        m_LocalElmtSoln.m_QpGradUold[i](1)+=dN(1)*uold;
                        m_LocalElmtSoln.m_QpGradUold[i](2)+=dN(2)*uold;
                        m_LocalElmtSoln.m_QpGradUold[i](3)+=dN(3)*uold;
                        //
                        m_LocalElmtSoln.m_QpGradU[i](1)+=dN(1)*u;
                        m_LocalElmtSoln.m_QpGradU[i](2)+=dN(2)*u;
                        m_LocalElmtSoln.m_QpGradU[i](3)+=dN(3)*u;

                        m_LocalElmtSoln.m_QpGradV[i](1)+=dN(1)*v;
                        m_LocalElmtSoln.m_QpGradV[i](2)+=dN(2)*v;
                        m_LocalElmtSoln.m_QpGradV[i](3)+=dN(3)*v;
                    }
                }// end-of-sub-element-dofs-loop


                //***********************************************************
                //*** for materials (UMAT) and elements/models (UEL)
                //***********************************************************
                t_MateSystem.runBulkMateLibs(t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_MateType,
                                             t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_JsonParams,
                                             m_LocalElmtInfo,
                                             m_LocalElmtSoln);

            }//end-of-sub-element-loop
            // execute the local projection
            localProjectionAction(m_BulkElmtNodesNum,m_ElmtConn,JxW,t_FE.m_BulkShp,t_MateSystem.m_MaterialContainer,Data);

        }// end-of-qpoints-loop
    }// end-of-element-loop

    // release the ghost array pointers before destroying the ghost copies
    t_SolnSystem.m_Uolder.restoreGhostArrayRead(&_aUolder);
    t_SolnSystem.m_Uold.restoreGhostArrayRead(&_aUold);
    t_SolnSystem.m_Ucurrent.restoreGhostArrayRead(&_aU);
    t_SolnSystem.m_V.restoreGhostArrayRead(&_aV);
    t_SolnSystem.m_A.restoreGhostArrayRead(&_aA);

    t_SolnSystem.m_Ucurrent.destroyGhostCopy();//
    t_SolnSystem.m_Uold.destroyGhostCopy();
    t_SolnSystem.m_Uolder.destroyGhostCopy();

    // for the velocity and acceleration
    t_SolnSystem.m_V.destroyGhostCopy();
    t_SolnSystem.m_A.destroyGhostCopy();

    //******************************************************
    //*** after the local projection(in the element loop)
    //*** one needs to execute the global projection behavior
    //******************************************************
    // execute the global projection action
    globalProjectionAction(t_FECell,Data);

}