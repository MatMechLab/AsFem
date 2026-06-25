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
//+++ Date   : 2022.08.22
//+++ Purpose: Write the result to the vtu (VTK unstructured grid) file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "OutputSystem/VTUWriter.h"
#include <charconv>
#include <string>

void VTUWriter::saveResults(const string &t_filename,
                            const FECell &t_fecell,
                            const DofHandler &t_dofHandler,
                            SolutionSystem &t_solution,
                            ProjectionSystem &t_projection){
    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);

    if(t_dofHandler.getActiveDofs()){}

    t_solution.m_Ucurrent.makeGhostCopy();

    for (auto &it:t_projection.getProjectionDataRef().m_ScalarProjMateVecList) it.makeGhostCopy();
    for (auto &it:t_projection.getProjectionDataRef().m_VectorProjMateVecList) it.makeGhostCopy();
    for (auto &it:t_projection.getProjectionDataRef().m_Rank2ProjMateVecList) it.makeGhostCopy();
    for (auto &it:t_projection.getProjectionDataRef().m_Rank4ProjMateVecList) it.makeGhostCopy();

    if(m_rank==0){
        std::ofstream out;
        out.open(t_filename,std::ios::out|std::ios::binary);
        if(!out.is_open()){
            MessagePrinter::printErrorTxt("can\'t create/open "+t_filename+", please make sure you have the write permission");
            MessagePrinter::exitAsFem();
        }

        const int nNodes=t_fecell.getFECellNodesNum();
        const int nElmts=t_fecell.getFECellBulkElmtsNum();

        // Fast formatted output: format each value with std::to_chars (no locale / no format-string
        // parsing, far faster than ostream<< number formatting). scientific-with-precision-6 is
        // byte-for-byte identical to the previous "%.6e"/scientific<<setprecision(6) output, and the
        // integer form is unique, so the resulting file is bit-identical. Values are accumulated into
        // an in-memory buffer and written out in ~1MB chunks (one bulk write instead of hundreds of
        // thousands of per-value stream insertions).
        std::string sbuf;
        sbuf.reserve(1u<<21);// 2 MB
        char nb[48];
        auto FLUSH=[&](){ if(!sbuf.empty()){ out.write(sbuf.data(),static_cast<std::streamsize>(sbuf.size())); sbuf.clear(); } };
        auto CHECK=[&](){ if(sbuf.size()>=(1u<<20)) FLUSH(); };// flush at ~1 MB
        auto D=[&](const double v){ auto r=std::to_chars(nb,nb+sizeof(nb),v,std::chars_format::scientific,6); sbuf.append(nb,static_cast<size_t>(r.ptr-nb)); };
        auto LI=[&](const long v){ auto r=std::to_chars(nb,nb+sizeof(nb),v); sbuf.append(nb,static_cast<size_t>(r.ptr-nb)); };

        sbuf+="<?xml version=\"1.0\"?>\n";
        sbuf+="<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n";
        sbuf+="<UnstructuredGrid>\n";
        sbuf+="<Piece NumberOfPoints=\""; LI(nNodes); sbuf+="\" NumberOfCells=\""; LI(nElmts); sbuf+="\">\n";
        sbuf+="<Points>\n";
        sbuf+="<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";

        //*****************************
        // print out node coordinates
        for (int i=1;i<=nNodes;i++){
            D(t_fecell.getFECellIthNodeJthCoord(i,1)); sbuf+=' ';
            D(t_fecell.getFECellIthNodeJthCoord(i,2)); sbuf+=' ';
            D(t_fecell.getFECellIthNodeJthCoord(i,3)); sbuf+='\n';
            CHECK();
        }
        sbuf+="</DataArray>\n";
        sbuf+="</Points>\n";

        //***************************************
        //*** For cell connectivity
        //***************************************
        sbuf+="<Cells>\n";
        sbuf+="<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (int e=1;e<=nElmts;e++){
            const int nn=t_fecell.getFECellIthBulkElmtNodesNum(e);
            for (int j=1;j<=nn;j++){
                LI(t_fecell.getFECellIthBulkElmtJthNodeID(e,j)-1); sbuf+=' ';
            }
            sbuf+='\n';
            CHECK();
        }
        sbuf+="</DataArray>\n";

        //***************************************
        //*** For offset
        //***************************************
        sbuf+="<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        long offset=0;
        for (int e=1;e<=nElmts;e++){
            offset+=t_fecell.getFECellIthBulkElmtNodesNum(e);
            LI(offset); sbuf+='\n';
            CHECK();
        }
        sbuf+="</DataArray>\n";

        //***************************************
        //*** For cell types
        //***************************************
        sbuf+="<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for (int e=1;e<=nElmts;e++){
            LI(t_fecell.getFECellIthBulkElmtVTKCellType(e)); sbuf+='\n';
            CHECK();
        }
        sbuf+="</DataArray>\n";
        sbuf+="</Cells>\n";

        //***************************************
        //*** For the PointData header
        //***************************************
        int i,j,k,iInd;
        string dofname;
        string ScalarName,VectorName,Rank2Name,Rank4Name,TensorName;

        ScalarName="<PointData Scalar=\"";
        for (j=1;j<=t_dofHandler.getMaxDofsPerNode();j++){
            ScalarName+=t_dofHandler.getIthDofName(j)+" ";
        }
        for(j=1;j<=t_projection.getScalarMaterialNum();j++){
            ScalarName+=t_projection.getIthScalarMateName(j)+" ";
        }
        ScalarName+="\" ";

        VectorName.clear();
        if(t_projection.getVectorMaterialNum()){
            VectorName="Vector=\"";
            for(j=1;j<=t_projection.getVectorMaterialNum();j++){
                VectorName+=t_projection.getIthVectorMateName(j)+" ";
            }
            VectorName+="\" ";
        }

        TensorName.clear();
        if(t_projection.getRank2MaterialNum()||t_projection.getRank4MaterialNum()){
            TensorName="Tensor=\"";
            for(j=1;j<=t_projection.getRank2MaterialNum();j++){
                TensorName+=t_projection.getIthRank2MateName(j)+" ";
            }
            for(j=1;j<=t_projection.getRank4MaterialNum();j++){
                TensorName+=t_projection.getIthRank4MateName(j)+" ";
            }
            TensorName+="\" ";
        }

        sbuf+=ScalarName; sbuf+=VectorName; sbuf+=TensorName; sbuf+=">\n";

        //**************************************
        //*** for solution (nodal dofs) output
        //**************************************
        const PetscScalar *usol=t_solution.m_Ucurrent.getGhostArrayRead();
        for (j=1;j<=t_dofHandler.getMaxDofsPerNode();j++){
            dofname=t_dofHandler.getIthDofName(j);
            sbuf+="<DataArray type=\"Float64\" Name=\""; sbuf+=dofname; sbuf+="\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (i=1;i<=nNodes;i++){
                iInd=t_dofHandler.getIthNodeJthDofID(i,j);
                D(usol[iInd-1]); sbuf+='\n';
                CHECK();
            }
            sbuf+="</DataArray>\n\n";
        }
        t_solution.m_Ucurrent.restoreGhostArrayRead(&usol);

        int nproj;
        //**************************************
        //*** for projected scalar output
        //**************************************
        nproj=t_projection.getScalarMaterialNum();
        for (j=1;j<=nproj;j++){
            const PetscScalar *parr=t_projection.getProjectionDataRef().m_ScalarProjMateVecList[j-1].getGhostArrayRead();
            dofname=t_projection.getIthScalarMateName(j);
            sbuf+="<DataArray type=\"Float64\" Name=\""; sbuf+=dofname; sbuf+="\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (i=1;i<=nNodes;i++){
                iInd=(i-1)*(1+1)+2;
                D(parr[iInd-1]); sbuf+='\n';
                CHECK();
            }
            sbuf+="</DataArray>\n\n";
            t_projection.getProjectionDataRef().m_ScalarProjMateVecList[j-1].restoreGhostArrayRead(&parr);
        }

        //**************************************
        //*** for projected vector output
        //**************************************
        nproj=t_projection.getVectorMaterialNum();
        for (j=1;j<=nproj;j++){
            const PetscScalar *parr=t_projection.getProjectionDataRef().m_VectorProjMateVecList[j-1].getGhostArrayRead();
            dofname=t_projection.getIthVectorMateName(j);
            sbuf+="<DataArray type=\"Float64\" Name=\""; sbuf+=dofname; sbuf+="\"  NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (i=1;i<=nNodes;i++){
                for(k=1;k<=3;k++){
                    iInd=(i-1)*(1+3)+k+1;
                    D(parr[iInd-1]); sbuf+=' ';
                }
                sbuf+='\n';
                CHECK();
            }
            sbuf+="</DataArray>\n\n";
            t_projection.getProjectionDataRef().m_VectorProjMateVecList[j-1].restoreGhostArrayRead(&parr);
        }

        //**************************************
        //*** for projected rank-2 output
        //**************************************
        nproj=t_projection.getRank2MaterialNum();
        for (j=1;j<=nproj;j++){
            const PetscScalar *parr=t_projection.getProjectionDataRef().m_Rank2ProjMateVecList[j-1].getGhostArrayRead();
            dofname=t_projection.getIthRank2MateName(j);
            sbuf+="<DataArray type=\"Float64\" Name=\""; sbuf+=dofname; sbuf+="\"  NumberOfComponents=\"9\" format=\"ascii\">\n";
            for (i=1;i<=nNodes;i++){
                for(k=1;k<=9;k++){
                    iInd=(i-1)*(1+9)+k+1;
                    D(parr[iInd-1]); sbuf+=' ';
                }
                sbuf+='\n';
                CHECK();
            }
            sbuf+="</DataArray>\n\n";
            t_projection.getProjectionDataRef().m_Rank2ProjMateVecList[j-1].restoreGhostArrayRead(&parr);
        }

        //**************************************
        //*** for projected rank-4 output
        //**************************************
        nproj=t_projection.getRank4MaterialNum();
        for (j=1;j<=nproj;j++){
            const PetscScalar *parr=t_projection.getProjectionDataRef().m_Rank4ProjMateVecList[j-1].getGhostArrayRead();
            dofname=t_projection.getIthRank4MateName(j);
            sbuf+="<DataArray type=\"Float64\" Name=\""; sbuf+=dofname; sbuf+="\"  NumberOfComponents=\"36\" format=\"ascii\">\n";
            for (i=1;i<=nNodes;i++){
                for(k=1;k<=36;k++){
                    iInd=(i-1)*(1+36)+k+1;
                    D(parr[iInd-1]); sbuf+=' ';
                }
                sbuf+='\n';
                CHECK();
            }
            sbuf+="</DataArray>\n\n";
            t_projection.getProjectionDataRef().m_Rank4ProjMateVecList[j-1].restoreGhostArrayRead(&parr);
        }

        //***************************************
        //*** End of output
        //***************************************
        sbuf+="</PointData>\n";
        sbuf+="</Piece>\n";
        sbuf+="</UnstructuredGrid>\n";
        sbuf+="</VTKFile>\n";
        FLUSH();

        out.close();
    }// end-of-master-rank-process

    t_solution.m_Ucurrent.destroyGhostCopy();
    for (auto &it:t_projection.getProjectionDataRef().m_ScalarProjMateVecList) it.destroyGhostCopy();
    for (auto &it:t_projection.getProjectionDataRef().m_VectorProjMateVecList) it.destroyGhostCopy();
    for (auto &it:t_projection.getProjectionDataRef().m_Rank2ProjMateVecList) it.destroyGhostCopy();
    for (auto &it:t_projection.getProjectionDataRef().m_Rank4ProjMateVecList) it.destroyGhostCopy();

}
