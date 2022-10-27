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
//+++ Date   : 2022.08.22
//+++ Purpose: Defines the abstract class for result output
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "OutputSystem/VTUWriter.h"

void VTUWriter::saveResults(const string &t_filename,
                            const Mesh &t_mesh,
                            const DofHandler &t_dofHandler,
                            SolutionSystem &t_solution,
                            ProjectionSystem &t_projection){
    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);

    t_solution.m_u_current.makeGhostCopy();
    t_projection.getProjectionDataRef().m_proj_scalarmate_vec.makeGhostCopy();
    t_projection.getProjectionDataRef().m_proj_vectormate_vec.makeGhostCopy();
    t_projection.getProjectionDataRef().m_proj_rank2mate_vec.makeGhostCopy();
    t_projection.getProjectionDataRef().m_proj_rank4mate_vec.makeGhostCopy();

    if(m_rank==0){
        std::ofstream out;
        out.open(t_filename,std::ios::out);
        if(!out.is_open()){
            MessagePrinter::printErrorTxt("can\'t create/open "+t_filename+", please make sure you have the write permission");
            MessagePrinter::exitAsFem();
        }

        int i,j,k,iInd,e;

        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n";
        out << "<UnstructuredGrid>\n";
        out << "<Piece NumberOfPoints=\"" << t_mesh.getBulkMeshNodesNum() << "\" NumberOfCells=\"" << t_mesh.getBulkMeshBulkElmtsNum() << "\">\n";
        out << "<Points>\n";
        out << "<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";

        //*****************************
        // print out node coordinates
        out <<std::scientific << std::setprecision(6);
        for (i = 1; i <= t_mesh.getBulkMeshNodesNum(); i++){
            out << t_mesh.getBulkMeshIthNodeJthCoord0(i, 1) << " ";
            out << t_mesh.getBulkMeshIthNodeJthCoord0(i, 2) << " ";
            out << t_mesh.getBulkMeshIthNodeJthCoord0(i, 3) << "\n";
        }
        out << "</DataArray>\n";
        out << "</Points>\n";

        //***************************************
        //*** For cell information
        //***************************************
        out << "<Cells>\n";
        out << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (e = 1; e <= t_mesh.getBulkMeshBulkElmtsNum(); e++){
            for (j = 1; j <= t_mesh.getBulkMeshIthBulkElmtNodesNum(e); j++){
                out << t_mesh.getBulkMeshIthBulkElmtJthNodeID(e,j) - 1 << " ";
            }
            out << "\n";
        }
        out << "</DataArray>\n";

        //***************************************
        //*** For offset
        //***************************************
        out << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int offset = 0;
        for (e = 1; e <= t_mesh.getBulkMeshBulkElmtsNum(); e++){
            offset += t_mesh.getBulkMeshIthBulkElmtNodesNum(e);
            out << offset << "\n";
        }
        out << "</DataArray>\n";

        //***************************************
        //*** For connectivity
        //***************************************
        out << "<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for (e = 1; e <= t_mesh.getBulkMeshBulkElmtsNum(); e++){
            out << t_mesh.getBulkMeshBulkElmtVTKCellType() << "\n";
        }
        out << "</DataArray>\n";
        out << "</Cells>\n";


        //***************************************
        //*** For solutions
        //***************************************
        string dofname;
        double value;
        string ScalarName,VectorName,Rank2Name,Rank4Name,TensorName;

        ScalarName="<PointData Scalar=\"";
        for (j = 1;j<=t_dofHandler.getMaxDofsPerNode();j++){
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

        out << ScalarName<< VectorName<< TensorName<<">\n";

        //**************************************
        //*** for solution output
        //**************************************
        for (j = 1;j<=t_dofHandler.getMaxDofsPerNode();j++){
            dofname =t_dofHandler.getIthDofName(j);
            out<<"<DataArray type=\"Float64\" Name=\"" << dofname << "\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            out<<std::scientific<<std::setprecision(6);
            for (i = 1; i <= t_mesh.getBulkMeshNodesNum(); i++){
                iInd = t_dofHandler.getIthNodeJthDofID(i,j);
                value=t_solution.m_u_current.getIthValueFromGhost(iInd);
                out << value << "\n";
            }
            out << "</DataArray>\n\n";
        }


        int nproj;
        //**************************************
        //*** for projected scalar output
        //**************************************
        nproj=t_projection.getScalarMaterialNum();
        for (j = 1;j<=nproj;j++){
            dofname = t_projection.getIthScalarMateName(j);
            out<<"<DataArray type=\"Float64\" Name=\"" << dofname << "\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            out<<std::scientific<<std::setprecision(6);
            for (i = 1; i <= t_mesh.getBulkMeshNodesNum(); i++){
                iInd = (i-1)*(1+nproj)+j+1;
                value=t_projection.getProjectionDataRef().m_proj_scalarmate_vec.getIthValueFromGhost(iInd);
                out << value << "\n";
            }
            out << "</DataArray>\n\n";
        }

        //**************************************
        //*** for projected vector output
        //**************************************
        nproj=t_projection.getVectorMaterialNum();
        for (j = 1;j<=nproj;j++){
            dofname = t_projection.getIthVectorMateName(j);
            out<<"<DataArray type=\"Float64\" Name=\"" << dofname << "\"  NumberOfComponents=\"3\" format=\"ascii\">\n";
            out<<std::scientific<<std::setprecision(6);
            for (i = 1; i <= t_mesh.getBulkMeshNodesNum(); i++){
                for(k=1;k<=3;k++){
                    iInd = (i-1)*(1+nproj*3)+3*(j-1)+k+1;
                    value=t_projection.getProjectionDataRef().m_proj_vectormate_vec.getIthValueFromGhost(iInd);
                    out << value << " ";
                }
                out <<"\n";
            }
            out << "</DataArray>\n\n";
        }

        //**************************************
        //*** for projected rank-2 output
        //**************************************
        nproj=t_projection.getRank2MaterialNum();
        for (j = 1;j<=nproj;j++){
            dofname = t_projection.getIthRank2MateName(j);
            out<<"<DataArray type=\"Float64\" Name=\"" << dofname << "\"  NumberOfComponents=\"9\" format=\"ascii\">\n";
            out<<std::scientific<<std::setprecision(6);
            for (i = 1; i <= t_mesh.getBulkMeshNodesNum(); i++){
                for(k=1;k<=9;k++){
                    iInd = (i-1)*(1+nproj*9)+9*(j-1)+k+1;
                    value=t_projection.getProjectionDataRef().m_proj_rank2mate_vec.getIthValueFromGhost(iInd);
                    out << value << " ";
                }
                out <<"\n";
            }
            out << "</DataArray>\n\n";
        }

        //**************************************
        //*** for projected rank-4 output
        //**************************************
        nproj=t_projection.getRank4MaterialNum();
        for (j = 1;j<=nproj;j++){
            dofname = t_projection.getIthRank4MateName(j);
            out<<"<DataArray type=\"Float64\" Name=\"" << dofname << "\"  NumberOfComponents=\"36\" format=\"ascii\">\n";
            out<<std::scientific<<std::setprecision(6);
            for (i = 1; i <= t_mesh.getBulkMeshNodesNum(); i++){
                for(k=1;k<=36;k++){
                    iInd = (i-1)*(1+nproj*36)+36*(j-1)+k+1;
                    value=t_projection.getProjectionDataRef().m_proj_rank4mate_vec.getIthValueFromGhost(iInd);
                    out << value << " ";
                }
                out <<"\n";
            }
            out << "</DataArray>\n\n";
        }

        

        //***************************************
        //*** End of output
        //***************************************
        out << "</PointData>\n";
        out << "</Piece>\n";
        out << "</UnstructuredGrid>\n";
        out << "</VTKFile>" << endl;

        out.close();

    }// end-of-master-rank-process

    t_solution.m_u_current.destroyGhostCopy();
    t_projection.getProjectionDataRef().m_proj_scalarmate_vec.destroyGhostCopy();
    t_projection.getProjectionDataRef().m_proj_vectormate_vec.destroyGhostCopy();
    t_projection.getProjectionDataRef().m_proj_rank2mate_vec.destroyGhostCopy();
    t_projection.getProjectionDataRef().m_proj_rank4mate_vec.destroyGhostCopy();

}