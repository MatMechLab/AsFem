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
//+++ Date   : 2025.01.24
//+++ Purpose: Implement full least square projection from integration
//+++          points to nodal points
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "ProjectionSystem/FullLeastSquareProjection.h"

FullLeastSquareProjection::FullLeastSquareProjection(){
    IsFirstSetup=true;
}
void FullLeastSquareProjection::initMyProjection(const FECell &t_FECell, const DofHandler &t_DofHandler) {
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
    IsFirstSetup=true;
}

//******************************************************
//*** for local projection action
//******************************************************
void FullLeastSquareProjection::localProjectionAction(const int &NodesNum,
                                                      const int &ElmtID,
                                                      const int &QPointsNum,
                                                      Eigen::MatrixXd &ALeft,
                                                      Eigen::MatrixXd &ARight,
                                                      const vector<int> &ElConn,
                                                      const double &Volume,
                                                      const SolutionSystem &SolnSystem,
                                                      ProjectionData &Data){

    if(Data.m_ScalarProjMateNum){
        projectLocalScalarMate2Global(NodesNum,ElmtID,QPointsNum,ALeft,ARight,ElConn,Volume,SolnSystem,Data.m_ScalarProjMateNum,Data.m_ScalarProjMateNameList,Data.m_ScalarProjMateVecList);
    }
    if(Data.m_VectorProjMateNum){
        projectLocalVectorMate2Global(NodesNum,ElmtID,QPointsNum,ALeft,ARight,ElConn,Volume,SolnSystem,Data.m_VectorProjMateNum,Data.m_VectorProjMateNamelist,Data.m_VectorProjMateVecList);
    }
    if(Data.m_Rank2ProjMateNum){
        projectLocalRank2Mate2Global(NodesNum,ElmtID,QPointsNum,ALeft,ARight,ElConn,Volume,SolnSystem,Data.m_Rank2ProjMateNum,Data.m_Rank2ProjMateNameList,Data.m_Rank2ProjMateVecList);
    }
    if(Data.m_Rank4ProjMateNum){
        projectLocalRank4Mate2Global(NodesNum,ElmtID,QPointsNum,ALeft,ARight,ElConn,Volume,SolnSystem,Data.m_Rank4ProjMateNum,Data.m_Rank4ProjMateNameList,Data.m_Rank4ProjMateVecList);
    }
}
//******************************************************
void FullLeastSquareProjection::projectLocalScalarMate2Global(const int &NodesNum,
                                                              const int &ElmtID,
                                                              const int &QPointsNum,
                                                              Eigen::MatrixXd &ALeft,
                                                              Eigen::MatrixXd &ARight,
                                                              const vector<int> &ElConn,
                                                              const double &Volume,
                                                              const SolutionSystem &SolnSystem,
                                                              const int &ProjNum,
                                                              const vector<string> &ScalarMateNameList,
                                                              vector<Vector> &ScalarProjMateVecList){
    if(ProjNum!=static_cast<int>(ScalarMateNameList.size())){
        MessagePrinter::printErrorTxt("the number of projected scalar material does not match with the given one, please check your code");
        MessagePrinter::exitAsFem();
    }

    bool HasName;
    int j,k,jInd,iInd,qp;
    for(k=1;k<=ProjNum;k++){
        HasName=false;
        for(const auto &it:SolnSystem.m_QpointsScalarMaterials_Local[(ElmtID-1)*QPointsNum+1-1]){
            if(it.first==ScalarMateNameList[k-1]){
                RHSVec.setZero();
                QpSolnVec.setZero();
                NodalSolnVec.setZero();
                for (qp=1;qp<=QPointsNum;qp++) {
                    QpSolnVec.coeffRef(qp-1)=SolnSystem.m_QpointsScalarMaterials_Local[(ElmtID-1)*QPointsNum+qp-1].at(it.first);
                }
                RHSVec=ARight*QpSolnVec;
                NodalSolnVec=ALeft.fullPivLu().solve(RHSVec);
                HasName=true;
                for(j=1;j<=NodesNum;j++) {
                    iInd=ElConn[j-1];
                    jInd=(iInd-1)*(1+1)+1;
                    ScalarProjMateVecList[k-1].addValue(jInd,Volume);
                    jInd=(iInd-1)*(1+1)+2;
                    ScalarProjMateVecList[k-1].addValue(jInd,NodalSolnVec.coeff(j-1)*Volume);
                }
                break;
            }
        }
        if(!HasName){
            MessagePrinter::printErrorTxt("can\'t find the scalar material ("+ScalarMateNameList[k-1]+") for projection, please check either your input file or your material code");
            MessagePrinter::exitAsFem();
        }
    }
}
//******************************************************
void FullLeastSquareProjection::projectLocalVectorMate2Global(const int &NodesNum,
                                                              const int &ElmtID,
                                                              const int &QPointsNum,
                                                              Eigen::MatrixXd &ALeft,
                                                              Eigen::MatrixXd &ARight,
                                                              const vector<int> &ElConn,
                                                              const double &Volume,
                                                              const SolutionSystem &SolnSystem,
                                                              const int &ProjNum,
                                                              const vector<string> &VectorMateNameList,
                                                              vector<Vector> &VectorProjMateVecList){
    if(ProjNum!=static_cast<int>(VectorMateNameList.size())){
        MessagePrinter::printErrorTxt("the number of projected vector material does not match with the given one, please check your code");
        MessagePrinter::exitAsFem();
    }
    bool HasName;
    int j,jInd,iInd;
    for(int k=1;k<=ProjNum;k++){
        HasName=false;
        for(const auto &it:SolnSystem.m_QpointsVectorMaterials_Local[(ElmtID-1)*QPointsNum+1-1]){
            if(it.first==VectorMateNameList[k-1]){
                HasName=true;
                for(j=1;j<=NodesNum;j++) {
                    iInd=ElConn[j-1];
                    jInd=(iInd-1)*(1+3)+1;
                    VectorProjMateVecList[k-1].addValue(jInd,Volume);
                }
                for (int component=1;component<=3;component++) {
                    RHSVec.setZero();
                    QpSolnVec.setZero();
                    NodalSolnVec.setZero();
                    for (int qp=1;qp<=QPointsNum;qp++) {
                        QpSolnVec.coeffRef(qp-1)=SolnSystem.m_QpointsVectorMaterials_Local[(ElmtID-1)*QPointsNum+qp-1].at(it.first)(component);
                    }
                    RHSVec=ARight*QpSolnVec;
                    NodalSolnVec=ALeft.fullPivLu().solve(RHSVec);
                    for (j=1;j<=NodesNum;j++) {
                        iInd=ElConn[j-1];
                        jInd=(iInd-1)*(1+3)+1+component;
                        VectorProjMateVecList[k-1].addValue(jInd,NodalSolnVec.coeff(j-1)*Volume);
                    }
                }
                break;
            }
        }
        if(!HasName){
            MessagePrinter::printErrorTxt("can\'t find the vector material ("+VectorMateNameList[k-1]+") for projection, please check either your input file or your material code");
            MessagePrinter::exitAsFem();
        }
    }
}
//******************************************************
void FullLeastSquareProjection::projectLocalRank2Mate2Global(const int &NodesNum,
                                                             const int &ElmtID,
                                                             const int &QPointsNum,
                                                             Eigen::MatrixXd &ALeft,
                                                             Eigen::MatrixXd &ARight,
                                                             const vector<int> &ElConn,
                                                             const double &Volume,
                                                             const SolutionSystem &SolnSystem,
                                                             const int &ProjNum,
                                                             const vector<string> &Rank2MateNameList,
                                                             vector<Vector> &Rank2ProjMateVecList){
    if(ProjNum!=static_cast<int>(Rank2MateNameList.size())){
        MessagePrinter::printErrorTxt("the number of projected rank-2 material does not match with the given one, please check your code");
        MessagePrinter::exitAsFem();
    }
    bool HasName;
    int j,k,jInd,iInd;
    int i1,j1,ii;
    for(k=1;k<=ProjNum;k++){
        HasName=false;
        for(const auto &it:SolnSystem.m_QpointsRank2Materials_Local[(ElmtID-1)*QPointsNum+1-1]){
            if(it.first==Rank2MateNameList[k-1]){
                HasName=true;
                for(j=1;j<=NodesNum;j++) {
                    iInd=ElConn[j-1];
                    jInd=(iInd-1)*(9+1)+0+1;
                    Rank2ProjMateVecList[k-1].addValue(jInd,Volume);
                }

                ii=0;
                for(i1=1;i1<=3;i1++){
                    for(j1=1;j1<=3;j1++){
                        RHSVec.setZero();
                        QpSolnVec.setZero();
                        NodalSolnVec.setZero();
                        for (int qp=1;qp<=QPointsNum;qp++) {
                            QpSolnVec.coeffRef(qp-1)=SolnSystem.m_QpointsRank2Materials_Local[(ElmtID-1)*QPointsNum+qp-1].at(it.first)(i1,j1);
                        }
                        RHSVec=ARight*QpSolnVec;
                        NodalSolnVec=ALeft.fullPivLu().solve(RHSVec);
                        ii+=1;
                        for(j=1;j<=NodesNum;j++) {
                            iInd=ElConn[j-1];
                            jInd=(iInd-1)*(9+1)+ii+1;
                            Rank2ProjMateVecList[k-1].addValue(jInd,NodalSolnVec.coeff(j-1)*Volume);
                        }
                    }
                }
                break;
            }
        }
        if(!HasName){
            MessagePrinter::printErrorTxt("can\'t find the rank-2 material ("+Rank2MateNameList[k-1]+") for projection, please check either your input file or your material code");
            MessagePrinter::exitAsFem();
        }
    }
}
//******************************************************
void FullLeastSquareProjection::projectLocalRank4Mate2Global(const int &NodesNum,
                                                             const int &ElmtID,
                                                             const int &QPointsNum,
                                                             Eigen::MatrixXd &ALeft,
                                                             Eigen::MatrixXd &ARight,
                                                             const vector<int> &ElConn,
                                                             const double &Volume,
                                                             const SolutionSystem &SolnSystem,
                                                             const int &ProjNum,
                                                             const vector<string> &Rank4MateNameList,
                                                             vector<Vector> &Rank4ProjMateVecList){
    if(ProjNum!=static_cast<int>(Rank4MateNameList.size())){
        MessagePrinter::printErrorTxt("the number of projected rank-4 material does not match with the given one, please check your code");
        MessagePrinter::exitAsFem();
    }
    bool HasName;
    int j,k,jInd,iInd;
    int i1,j1,ii;
    for(k=1;k<=ProjNum;k++){
        HasName=false;
        for(const auto &it:SolnSystem.m_QpointsRank4Materials_Local[(ElmtID-1)*QPointsNum+1-1]){
            if(it.first==Rank4MateNameList[k-1]){
                HasName=true;
                ii=0;
                for(j=1;j<=NodesNum;j++) {
                    iInd=ElConn[j-1];
                    jInd=(iInd-1)*(36+1)+ii+1;
                    Rank4ProjMateVecList[k-1].addValue(jInd,Volume);
                }
                for(i1=1;i1<=6;i1++){
                    for(j1=1;j1<=6;j1++){
                        RHSVec.setZero();
                        QpSolnVec.setZero();
                        NodalSolnVec.setZero();
                        for (int qp=1;qp<=QPointsNum;qp++) {
                            QpSolnVec.coeffRef(qp-1)=SolnSystem.m_QpointsRank4Materials_Local[(ElmtID-1)*QPointsNum+qp-1].at(it.first).getVoigtComponent(i1,j1);
                        }
                        RHSVec=ARight*QpSolnVec;
                        NodalSolnVec=ALeft.fullPivLu().solve(RHSVec);

                        ii+=1;
                        for(j=1;j<=NodesNum;j++) {
                            iInd=ElConn[j-1];
                            jInd=(iInd-1)*(36+1)+ii+1;
                            Rank4ProjMateVecList[k-1].addValue(jInd,NodalSolnVec.coeff(j-1)*Volume);
                        }
                    }
                }
                break;
            }
        }
        if(!HasName){
            MessagePrinter::printErrorTxt("can\'t find the rank-4 material ("+Rank4MateNameList[k-1]+") for projection, please check either your input file or your material code");
            MessagePrinter::exitAsFem();
        }
    }
}
//******************************************************
//*** for global projection action
//******************************************************
void FullLeastSquareProjection::globalProjectionAction(const FECell &t_FECell,
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

    for (const auto &GlobalNodeID:t_FECell.getLocalFECellNodeIDsCopy()) {
        //****************************************
        // for scalar materials
        //****************************************
        nproj = Data.m_ScalarProjMateNum;
        if (nproj) {
            for (j = 1; j <= nproj; j++) {
                iInd = (GlobalNodeID - 1) * (1 + 1) + 1;
                weight = Data.m_ScalarProjMateVecList[j-1].getIthValueFromGhost(iInd);
                iInd = (GlobalNodeID - 1) * (1 + 1) + 2;
                value = Data.m_ScalarProjMateVecList[j-1].getIthValueFromGhost(iInd);
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
                weight = Data.m_VectorProjMateVecList[j-1].getIthValueFromGhost(iInd);
                for (k = 1; k <= 3; k++) {
                    iInd = (GlobalNodeID - 1) * (3 + 1)+ k + 1;
                    value = Data.m_VectorProjMateVecList[j-1].getIthValueFromGhost(iInd);
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
                weight = Data.m_Rank2ProjMateVecList[j-1].getIthValueFromGhost(iInd);
                for (k = 1; k <= 9; k++) {
                    iInd = (GlobalNodeID - 1) * (9 + 1) + k + 1;
                    value = Data.m_Rank2ProjMateVecList[j-1].getIthValueFromGhost(iInd);
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
                weight = Data.m_Rank2ProjMateVecList[j-1].getIthValueFromGhost(iInd);
                for (k = 1; k <= 36; k++) {
                    iInd = (GlobalNodeID - 1) * (36 + 1)+ k + 1;
                    value = Data.m_Rank2ProjMateVecList[j-1].getIthValueFromGhost(iInd);
                    if (abs(weight) > 1.0e-15) {
                        newvalue = value / weight;
                    } else {
                        newvalue = value /1.0e-8;
                    }
                    Data.m_Rank2ProjMateVecList[j-1].addValue(iInd, newvalue);
                }
            }
        }
    } //end-of-local-element-nodes-loop


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
void FullLeastSquareProjection::executeMyProjection(const FECell &t_FECell,
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
    double xi,eta,zeta,w,J,JxW,Volume;

    int SubElmtBlockID;
    int GlobalDofID,GlobalNodeID;

    m_LocalElmtInfo.m_Dim=t_FECell.getFECellMaxDim();
    m_LocalElmtInfo.m_NodesNum=t_FECell.getFECellNodesNumPerBulkElmt();

    vector<SingleMeshCell> MyLocalCellVec;

    MyLocalCellVec=t_FECell.getLocalBulkFECellVecCopy();
    m_BulkElmtNodesNum=t_FECell.getFECellNodesNumPerBulkElmt();

    QpointsNum=t_FE.m_BulkQpoints.getQPointsNum();
    m_LocalElmtInfo.m_QpointsNum=QpointsNum;

    if (IsFirstSetup) {
        A.resize(QpointsNum,m_BulkElmtNodesNum);
        AT.resize(m_BulkElmtNodesNum,QpointsNum);
        W.resize(QpointsNum,QpointsNum);

        Al.resize(m_BulkElmtNodesNum,m_BulkElmtNodesNum);
        Ar.resize(m_BulkElmtNodesNum,QpointsNum);

        QpSolnVec.resize(QpointsNum);
        NodalSolnVec.resize(m_BulkElmtNodesNum);
        RHSVec.resize(m_BulkElmtNodesNum);
        IsFirstSetup=false;
    }

    for (int e=1;e<=t_FECell.getLocalFECellBulkElmtsNum();e++) {
        m_LocalElmtInfo.m_Dim=MyLocalCellVec[e-1].Dim;
        m_LocalElmtInfo.m_NodesNum=MyLocalCellVec[e-1].NodesNumPerElmt;
        m_LocalElmtInfo.m_Dt=t_FECtrlInfo.Dt;
        m_LocalElmtInfo.m_T=t_FECtrlInfo.T;

        m_Nodes=MyLocalCellVec[e-1].ElmtNodeCoords;
        m_Nodes0=MyLocalCellVec[e-1].ElmtNodeCoords0;
        m_ElmtConn=MyLocalCellVec[e-1].ElmtConn;

        A.setZero();
        AT.setZero();
        W.setZero();
        Al.setZero();
        Ar.setZero();

        Volume=0.0;

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
            Volume+=1.0*JxW;

            W.coeffRef(qp-1,qp-1)=w;

            m_LocalElmtInfo.m_QpointID=qp;
            m_LocalElmtInfo.m_QpCoords0=0.0;
            for (int i=1;i<=m_LocalElmtInfo.m_NodesNum;i++) {
                m_LocalElmtInfo.m_QpCoords0(1)+=t_FE.m_BulkShp.shape_value(i)*m_Nodes0(i,1);
                m_LocalElmtInfo.m_QpCoords0(2)+=t_FE.m_BulkShp.shape_value(i)*m_Nodes0(i,2);
                m_LocalElmtInfo.m_QpCoords0(3)+=t_FE.m_BulkShp.shape_value(i)*m_Nodes0(i,3);

                A.coeffRef(qp-1,i-1)=t_FE.m_BulkShp.shape_value(i);
                AT.coeffRef(i-1,qp-1)=t_FE.m_BulkShp.shape_value(i);
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
                        GlobalNodeID=MyLocalCellVec[e-1].ElmtConn[j-1];
                        GlobalDofID=t_DofHandler.getIthNodeJthDofID(GlobalNodeID,m_SubElmtDofIDs[i-1]);

                        m_LocalElmtSoln.m_QpUolder[i]+=t_FE.m_BulkShp.shape_value(j)*t_SolnSystem.m_Uolder.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpUold[i]  +=t_FE.m_BulkShp.shape_value(j)*t_SolnSystem.m_Uold.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpU[i]     +=t_FE.m_BulkShp.shape_value(j)*t_SolnSystem.m_Ucurrent.getIthValueFromGhost(GlobalDofID);

                        m_LocalElmtSoln.m_QpV[i]+=t_FE.m_BulkShp.shape_value(j)*t_SolnSystem.m_V.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpA[i]+=t_FE.m_BulkShp.shape_value(j)*t_SolnSystem.m_A.getIthValueFromGhost(GlobalDofID);

                        m_LocalElmtSoln.m_QpGradUolder[i](1)+=t_FE.m_BulkShp.shape_grad(j)(1)*t_SolnSystem.m_Uolder.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradUolder[i](2)+=t_FE.m_BulkShp.shape_grad(j)(2)*t_SolnSystem.m_Uolder.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradUolder[i](3)+=t_FE.m_BulkShp.shape_grad(j)(3)*t_SolnSystem.m_Uolder.getIthValueFromGhost(GlobalDofID);
                        //
                        m_LocalElmtSoln.m_QpGradUold[i](1)+=t_FE.m_BulkShp.shape_grad(j)(1)*t_SolnSystem.m_Uold.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradUold[i](2)+=t_FE.m_BulkShp.shape_grad(j)(2)*t_SolnSystem.m_Uold.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradUold[i](3)+=t_FE.m_BulkShp.shape_grad(j)(3)*t_SolnSystem.m_Uold.getIthValueFromGhost(GlobalDofID);
                        //
                        m_LocalElmtSoln.m_QpGradU[i](1)+=t_FE.m_BulkShp.shape_grad(j)(1)*t_SolnSystem.m_Ucurrent.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradU[i](2)+=t_FE.m_BulkShp.shape_grad(j)(2)*t_SolnSystem.m_Ucurrent.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradU[i](3)+=t_FE.m_BulkShp.shape_grad(j)(3)*t_SolnSystem.m_Ucurrent.getIthValueFromGhost(GlobalDofID);

                        m_LocalElmtSoln.m_QpGradV[i](1)+=t_FE.m_BulkShp.shape_grad(j)(1)*t_SolnSystem.m_V.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradV[i](2)+=t_FE.m_BulkShp.shape_grad(j)(2)*t_SolnSystem.m_V.getIthValueFromGhost(GlobalDofID);
                        m_LocalElmtSoln.m_QpGradV[i](3)+=t_FE.m_BulkShp.shape_grad(j)(3)*t_SolnSystem.m_V.getIthValueFromGhost(GlobalDofID);
                    }
                }// end-of-sub-element-dofs-loop


                //***********************************************************
                //*** for materials (UMAT) and elements/models (UEL)
                //***********************************************************
                t_MateSystem.runBulkMateLibs(t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_MateType,
                                             t_ElmtSystem.getIthBulkElmtBlock(SubElmtBlockID).m_JsonParams,
                                             m_LocalElmtInfo,
                                             m_LocalElmtSoln);

                t_SolnSystem.m_QpointsScalarMaterials_Local[(e-1)*QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getScalarMaterialsCopy();
                t_SolnSystem.m_QpointsVectorMaterials_Local[(e-1)*QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getVectorMaterialsCopy();
                t_SolnSystem.m_QpointsRank2Materials_Local[(e-1)*QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getRank2MaterialsCopy();
                t_SolnSystem.m_QpointsRank4Materials_Local[(e-1)*QpointsNum+qp-1]=t_MateSystem.m_MaterialContainer.getRank4MaterialsCopy();

            }//end-of-sub-element-loop
        }// end-of-qpoints-loop

        // execute the local projection, one must call this outside the qpoint loop
        // which means, this is the elemental operation !!!
        Al=AT*W*A;
        Ar=AT*W;
        if (Al.fullPivLu().rcond()<1.0e-12) {
            MessagePrinter::printErrorTxt("The extrapolation matrix is singular for the FullLeastSquare projection, please check your code");
            MessagePrinter::exitAsFem();
        }
        localProjectionAction(m_BulkElmtNodesNum,e,QpointsNum,Al,Ar,m_ElmtConn,Volume,t_SolnSystem,Data);
    }// end-of-element-loop

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