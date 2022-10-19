//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.26
//+++ Purpose: here we apply dirichlet boundary condition via the 
//+++          penalty method
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

void BCSystem::applyDirichletBC(const FECalcType &calctype,
                                const double &bcvalue,
                                const BCType &bctype,
                                const nlohmann::json &json,
                                const vector<int> &dofids,
                                const vector<string> &bcnamelist,
                                const Mesh &mesh,const DofHandler &dofhandler,
                                Vector &U,Vector &Ucopy,Vector &Uold,Vector &Uolder,
                                Vector &V,
                                SparseMatrix &AMATRIX,
                                Vector &RHS){
    //************************************
    //*** get rid of unused warnings 
    //************************************
    int i,j,k,e,iInd,nElmts;
    int rankne,eStart,eEnd;
    vector<int> globaldofids;
    globaldofids.resize(dofids.size(),0);

    m_local_elmtinfo.m_dofsnum=static_cast<int>(dofids.size());
    if(Uold.getSize()||Uolder.getSize()||V.getSize()) {}

    Ucopy.makeGhostCopy();
    Uold.makeGhostCopy();
    Uolder.makeGhostCopy();
    V.makeGhostCopy();

    MPI_Comm_size(PETSC_COMM_WORLD,&m_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);

    for(const auto &name:bcnamelist){
        nElmts=mesh.getBulkMeshElmtsNumViaPhyName(name);
        rankne=nElmts/m_size;
        eStart=m_rank*rankne;
        eEnd=(m_rank+1)*rankne;
        if(m_rank==m_size-1) eEnd=nElmts;
        m_local_elmtinfo.m_dim=mesh.getBulkMeshElmtDimViaPhyName(name);

        for(e=eStart;e<eEnd;e++){
            m_local_elmtinfo.m_nodesnum=mesh.getBulkMeshIthElmtNodesNumViaPhyName(name,e+1);

            for(i=1;i<=m_local_elmtinfo.m_nodesnum;i++){
                j=mesh.getBulkMeshIthElmtJthNodeIDViaPhyName(name,e+1,i);
                m_local_elmtinfo.m_gpCoords0(1)=mesh.getBulkMeshIthNodeJthCoord0(j,1);
                m_local_elmtinfo.m_gpCoords0(2)=mesh.getBulkMeshIthNodeJthCoord0(j,2);
                m_local_elmtinfo.m_gpCoords0(3)=mesh.getBulkMeshIthNodeJthCoord0(j,3);

                for(k=1;k<=m_local_elmtinfo.m_dofsnum;k++){
                    iInd=dofhandler.getIthNodeJthDofID(j,dofids[k-1]);
                    globaldofids[k-1]=iInd;

                    m_local_elmtsoln.m_gpU[k]=Ucopy.getIthValueFromGhost(iInd);
                    m_local_elmtsoln.m_gpUold[k]=Uold.getIthValueFromGhost(iInd);
                    m_local_elmtsoln.m_gpUolder[k]=Uolder.getIthValueFromGhost(iInd);
                    m_local_elmtsoln.m_gpV[k]=V.getIthValueFromGhost(iInd);
                }

                switch (bctype)
                {
                case BCType::DIRICHLETBC:
                    DirichletBC::computeBCValue(calctype,m_dirichlet_penalty,bcvalue,json,
                                                m_local_elmtinfo,m_local_elmtsoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                case BCType::USER1DIRICHLETBC:
                    User1DirichletBC::computeBCValue(calctype,m_dirichlet_penalty,bcvalue,json,
                                                m_local_elmtinfo,m_local_elmtsoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                case BCType::USER2DIRICHLETBC:
                    User2DirichletBC::computeBCValue(calctype,m_dirichlet_penalty,bcvalue,json,
                                                m_local_elmtinfo,m_local_elmtsoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                case BCType::USER3DIRICHLETBC:
                    User3DirichletBC::computeBCValue(calctype,m_dirichlet_penalty,bcvalue,json,
                                                m_local_elmtinfo,m_local_elmtsoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                case BCType::USER4DIRICHLETBC:
                    User4DirichletBC::computeBCValue(calctype,m_dirichlet_penalty,bcvalue,json,
                                                m_local_elmtinfo,m_local_elmtsoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                case BCType::USER5DIRICHLETBC:
                    User5DirichletBC::computeBCValue(calctype,m_dirichlet_penalty,bcvalue,json,
                                                m_local_elmtinfo,m_local_elmtsoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                case BCType::POISSON2DBENCHMARKBC:
                    Poisson2DBenchmarkBC::computeBCValue(calctype,m_dirichlet_penalty,bcvalue,json,
                                                m_local_elmtinfo,m_local_elmtsoln,globaldofids,
                                                U,AMATRIX,RHS);
                    break;
                default:
                    MessagePrinter::printErrorTxt("unsupported dirichlet type boundary condition in ApplyDirichletBC.cpp, plese check your input file or your code");
                    MessagePrinter::exitAsFem();
                    break;
                }
            }
        }
    }

    Ucopy.destroyGhostCopy();
    Uold.destroyGhostCopy();
    Uolder.destroyGhostCopy();
    V.destroyGhostCopy();

    // do the assemble
    U.assemble();
    if(calctype==FECalcType::COMPUTERESIDUAL) RHS.assemble();
    if(calctype==FECalcType::COMPUTEJACOBIAN) AMATRIX.assemble();

}