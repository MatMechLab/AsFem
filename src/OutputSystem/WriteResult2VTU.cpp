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
//+++ Date   : 2020.12.17
//+++ Purpose: write results to vtu file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "OutputSystem/OutputSystem.h"

void OutputSystem::WriteResult2VTU(const Mesh &mesh,const DofHandler &dofHandler,const Vec &U){
    MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);

    VecScatterCreateToAll(U,&_scatterU,&_Useq);
    VecScatterBegin(_scatterU,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterU,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);

    if(_rank == 0){
        _OutputFileName=_InputFileName.substr(0,_InputFileName.size()-2);// remove ".i" extension name
        _VTUFileName = _OutputFileName + ".vtu";
        ofstream _VTUFile;
        _VTUFile.open(_VTUFileName, ios::out);
        if (!_VTUFile.is_open()){
            string str="can\'t create a new vtu file(="+_VTUFileName+")!, please make sure you have write permission";
            MessagePrinter::PrintErrorTxt(str);
            MessagePrinter::AsFem_Exit();
        }
        int i,j,iInd,e;

        _OutputFileName=_VTUFileName;
        //****************************************
        //*** print out header information
        //****************************************
        _VTUFile << "<?xml version=\"1.0\"?>\n";
        _VTUFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
        _VTUFile << "<UnstructuredGrid>\n";
        _VTUFile << "<Piece NumberOfPoints=\"" << mesh.GetBulkMeshNodesNum() << "\" NumberOfCells=\"" << mesh.GetBulkMeshBulkElmtsNum() << "\">\n";
        _VTUFile << "<Points>\n";
        _VTUFile << "<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";

        //*****************************
        // print out node coordinates
        _VTUFile<<scientific << setprecision(6);
        for (i = 1; i <= mesh.GetBulkMeshNodesNum(); ++i){
            _VTUFile << mesh.GetBulkMeshIthNodeJthCoord(i, 1) << " ";
            _VTUFile << mesh.GetBulkMeshIthNodeJthCoord(i, 2) << " ";
            _VTUFile << mesh.GetBulkMeshIthNodeJthCoord(i, 3) << "\n";
        }
        _VTUFile << "</DataArray>\n";
        _VTUFile << "</Points>\n";

        //***************************************
        //*** For cell information
        //***************************************
        _VTUFile << "<Cells>\n";
        _VTUFile << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (e = 1; e <= mesh.GetBulkMeshBulkElmtsNum(); ++e){
            for (j = 1; j <= mesh.GetBulkMeshIthBulkElmtNodesNum(e); ++j){
                _VTUFile << mesh.GetBulkMeshIthBulkElmtJthNodeID(e, j) - 1 << " ";
            }
            _VTUFile << "\n";
        }
        _VTUFile << "</DataArray>\n";

        // for offset
        _VTUFile << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int offset = 0;
        for (e = 1; e <= mesh.GetBulkMeshBulkElmtsNum(); ++e){
            offset += mesh.GetBulkMeshIthBulkElmtNodesNum(e);
            _VTUFile << offset << "\n";
        }
        _VTUFile << "</DataArray>\n";

        // For connectivity
        _VTUFile << "<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for (e = 1; e <= mesh.GetBulkMeshBulkElmtsNum(); ++e){
            _VTUFile << mesh.GetBulkMeshIthBulkElmtVTKCellType(e) << "\n";
        }
        _VTUFile << "</DataArray>\n";
        _VTUFile << "</Cells>\n";

        // for point data
        _VTUFile << "<PointData Scalar=\"sol\">\n";

        string dofname;
        PetscScalar value;

        for (j = 1; j <= dofHandler.GetDofsNumPerNode(); ++j){
            dofname = dofHandler.GetIthDofName(j);
            _VTUFile << "<DataArray type=\"Float64\" Name=\"" << dofname << "\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            _VTUFile << scientific << setprecision(6);
            for (i = 1; i <= mesh.GetBulkMeshNodesNum(); ++i){
                iInd = dofHandler.GetIthNodeJthDofIndex(i, j) - 1;
                VecGetValues(_Useq,1,&iInd,&value);
                _VTUFile << value << "\n";
            }
            _VTUFile << "</DataArray>\n\n";
        }
        _VTUFile << "</PointData>\n";
        _VTUFile << "</Piece>\n";
        _VTUFile << "</UnstructuredGrid>\n";
        _VTUFile << "</VTKFile>" << endl;

        _VTUFile.close();
    }

    VecScatterDestroy(&_scatterU);
    VecDestroy(&_Useq);
}
void OutputSystem::WriteResult2VTU(const int &step, const Mesh &mesh, const DofHandler &dofHandler, const Vec &U){
    MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);

    VecScatterCreateToAll(U,&_scatterU,&_Useq);
    VecScatterBegin(_scatterU,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterU,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);

    if(_rank == 0){
        ostringstream ss;
        ss<<setfill('0')<<setw(8)<<step;
        _OutputFileName=_InputFileName.substr(0,_InputFileName.size()-2);// remove ".i" extension name
        _VTUFileName = _OutputFileName+"-"+ss.str()+ ".vtu";
        ofstream _VTUFile;
        _VTUFile.open(_VTUFileName, ios::out);
        if (!_VTUFile.is_open()){
            string str="can\'t create a new vtu file(="+_VTUFileName+")!, please make sure you have write permission";
            MessagePrinter::PrintErrorTxt(str);
            MessagePrinter::AsFem_Exit();
        }
        int i,j,iInd,e;

        _OutputFileName=_VTUFileName;
        //****************************************
        //*** print out header information
        //****************************************
        _VTUFile << "<?xml version=\"1.0\"?>\n";
        _VTUFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
        _VTUFile << "<UnstructuredGrid>\n";
        _VTUFile << "<Piece NumberOfPoints=\"" << mesh.GetBulkMeshNodesNum() << "\" NumberOfCells=\"" << mesh.GetBulkMeshBulkElmtsNum() << "\">\n";
        _VTUFile << "<Points>\n";
        _VTUFile << "<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";

        //*****************************
        // print out node coordinates
        _VTUFile<<scientific << setprecision(6);
        for (i = 1; i <= mesh.GetBulkMeshNodesNum(); ++i){
            _VTUFile << mesh.GetBulkMeshIthNodeJthCoord(i, 1) << " ";
            _VTUFile << mesh.GetBulkMeshIthNodeJthCoord(i, 2) << " ";
            _VTUFile << mesh.GetBulkMeshIthNodeJthCoord(i, 3) << "\n";
        }
        _VTUFile << "</DataArray>\n";
        _VTUFile << "</Points>\n";

        //***************************************
        //*** For cell information
        //***************************************
        _VTUFile << "<Cells>\n";
        _VTUFile << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (e = 1; e <= mesh.GetBulkMeshBulkElmtsNum(); ++e){
            for (j = 1; j <= mesh.GetBulkMeshIthBulkElmtNodesNum(e); ++j){
                _VTUFile << mesh.GetBulkMeshIthBulkElmtJthNodeID(e, j) - 1 << " ";
            }
            _VTUFile << "\n";
        }
        _VTUFile << "</DataArray>\n";

        // for offset
        _VTUFile << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int offset = 0;
        for (e = 1; e <= mesh.GetBulkMeshBulkElmtsNum(); ++e){
            offset += mesh.GetBulkMeshIthBulkElmtNodesNum(e);
            _VTUFile << offset << "\n";
        }
        _VTUFile << "</DataArray>\n";

        // For connectivity
        _VTUFile << "<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for (e = 1; e <= mesh.GetBulkMeshBulkElmtsNum(); ++e){
            _VTUFile << mesh.GetBulkMeshIthBulkElmtVTKCellType(e) << "\n";
        }
        _VTUFile << "</DataArray>\n";
        _VTUFile << "</Cells>\n";

        // for point data
        _VTUFile << "<PointData Scalar=\"sol\">\n";

        string dofname;
        PetscScalar value;

        for (j = 1; j <= dofHandler.GetDofsNumPerNode(); ++j){
            dofname = dofHandler.GetIthDofName(j);
            _VTUFile << "<DataArray type=\"Float64\" Name=\"" << dofname << "\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            _VTUFile << scientific << setprecision(6);
            for (i = 1; i <= mesh.GetBulkMeshNodesNum(); ++i){
                iInd = dofHandler.GetIthNodeJthDofIndex(i, j) - 1;
                VecGetValues(_Useq,1,&iInd,&value);
                _VTUFile << value << "\n";
            }
            _VTUFile << "</DataArray>\n\n";
        }
        _VTUFile << "</PointData>\n";
        _VTUFile << "</Piece>\n";
        _VTUFile << "</UnstructuredGrid>\n";
        _VTUFile << "</VTKFile>" << endl;

        _VTUFile.close();
    }

    VecScatterDestroy(&_scatterU);
    VecDestroy(&_Useq);
}
//*********************************************************************
void OutputSystem::WriteResult2VTU(const Mesh &mesh,const DofHandler &dofHandler,const Vec &U,const int &nProj,
                                   const vector<string> &projname,const Vec &Proj){
    MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);

    VecScatterCreateToAll(U,&_scatterU,&_Useq);
    VecScatterBegin(_scatterU,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterU,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterCreateToAll(Proj,&_scatterProj,&_ProjSeq);
    VecScatterBegin(_scatterProj,Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterProj,Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);

    if(_rank == 0){
        _OutputFileName=_InputFileName.substr(0,_InputFileName.size()-2);// remove ".i" extension name
        _VTUFileName = _OutputFileName + ".vtu";
        ofstream _VTUFile;
        _VTUFile.open(_VTUFileName, ios::out);
        if (!_VTUFile.is_open()){
            string str="can\'t create a new vtu file(="+_VTUFileName+")!, please make sure you have write permission";
            MessagePrinter::PrintErrorTxt(str);
            MessagePrinter::AsFem_Exit();
        }
        int i,j,iInd,e;

        _OutputFileName=_VTUFileName;

        //****************************************
        //*** print out header information
        //****************************************
        _VTUFile << "<?xml version=\"1.0\"?>\n";
        _VTUFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
        _VTUFile << "<UnstructuredGrid>\n";
        _VTUFile << "<Piece NumberOfPoints=\"" << mesh.GetBulkMeshNodesNum() << "\" NumberOfCells=\"" << mesh.GetBulkMeshBulkElmtsNum() << "\">\n";
        _VTUFile << "<Points>\n";
        _VTUFile << "<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";

        //*****************************
        // print out node coordinates
        _VTUFile<<scientific << setprecision(6);
        for (i = 1; i <= mesh.GetBulkMeshNodesNum(); ++i){
            _VTUFile << mesh.GetBulkMeshIthNodeJthCoord(i, 1) << " ";
            _VTUFile << mesh.GetBulkMeshIthNodeJthCoord(i, 2) << " ";
            _VTUFile << mesh.GetBulkMeshIthNodeJthCoord(i, 3) << "\n";
        }
        _VTUFile << "</DataArray>\n";
        _VTUFile << "</Points>\n";

        //***************************************
        //*** For cell information
        //***************************************
        _VTUFile << "<Cells>\n";
        _VTUFile << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (e = 1; e <= mesh.GetBulkMeshBulkElmtsNum(); ++e){
            for (j = 1; j <= mesh.GetBulkMeshIthBulkElmtNodesNum(e); ++j){
                _VTUFile << mesh.GetBulkMeshIthBulkElmtJthNodeID(e, j) - 1 << " ";
            }
            _VTUFile << "\n";
        }
        _VTUFile << "</DataArray>\n";

        // for offset
        _VTUFile << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int offset = 0;
        for (e = 1; e <= mesh.GetBulkMeshBulkElmtsNum(); ++e){
            offset += mesh.GetBulkMeshIthBulkElmtNodesNum(e);
            _VTUFile << offset << "\n";
        }
        _VTUFile << "</DataArray>\n";

        // For connectivity
        _VTUFile << "<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for (e = 1; e <= mesh.GetBulkMeshBulkElmtsNum(); ++e){
            _VTUFile << mesh.GetBulkMeshIthBulkElmtVTKCellType(e) << "\n";
        }
        _VTUFile << "</DataArray>\n";
        _VTUFile << "</Cells>\n";

        // for point data
        _VTUFile << "<PointData Scalar=\"sol\">\n";

        string dofname;
        PetscScalar value;

        for (j = 1; j <= dofHandler.GetDofsNumPerNode(); ++j){
            dofname = dofHandler.GetIthDofName(j);
            _VTUFile << "<DataArray type=\"Float64\" Name=\"" << dofname << "\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            _VTUFile << scientific << setprecision(6);
            for (i = 1; i <= mesh.GetBulkMeshNodesNum(); ++i){
                iInd = dofHandler.GetIthNodeJthDofIndex(i, j) - 1;
                VecGetValues(_Useq,1,&iInd,&value);
                _VTUFile << value << "\n";
            }
            _VTUFile << "</DataArray>\n\n";
        }

        //************************************
        //*** for projection value
        //************************************
        for(j=1;j<=nProj;++j){
            _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<projname[j-1]<<"\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            for(i=1;i<=mesh.GetBulkMeshNodesNum();++i){
                iInd =(i-1)*(1+nProj)+j;
                VecGetValues(_ProjSeq,1,&iInd,&value);
                _VTUFile<<scientific<<setprecision(6)<<value<<"\n";
            }
            _VTUFile<<"</DataArray>\n\n";
        }

        _VTUFile << "</PointData>\n";
        _VTUFile << "</Piece>\n";
        _VTUFile << "</UnstructuredGrid>\n";
        _VTUFile << "</VTKFile>" << endl;

        _VTUFile.close();
    }

    VecScatterDestroy(&_scatterU);
    VecDestroy(&_Useq);

    VecScatterDestroy(&_scatterProj);
    VecDestroy(&_ProjSeq);
}
//*******************************************************
void OutputSystem::WriteResult2VTU(const int &step, const Mesh &mesh, const DofHandler &dofHandler, const Vec &U,
                                   const int &nProj, const vector<string> &projname, const Vec &Proj){
    MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);

    VecScatterCreateToAll(U,&_scatterU,&_Useq);
    VecScatterBegin(_scatterU,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterU,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterCreateToAll(Proj,&_scatterProj,&_ProjSeq);
    VecScatterBegin(_scatterProj,Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterProj,Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);

    if(_rank == 0){
        ostringstream ss;
        ss<<setfill('0')<<setw(8)<<step;
        _OutputFileName=_InputFileName.substr(0,_InputFileName.size()-2);// remove ".i" extension name
        _VTUFileName = _OutputFileName+"-"+ss.str()+ ".vtu";
        ofstream _VTUFile;
        _VTUFile.open(_VTUFileName, ios::out);
        if (!_VTUFile.is_open()){
            string str="can\'t create a new vtu file(="+_VTUFileName+")!, please make sure you have write permission";
            MessagePrinter::PrintErrorTxt(str);
            MessagePrinter::AsFem_Exit();
        }
        int i,j,iInd,e;

        _OutputFileName=_VTUFileName;

        //****************************************
        //*** print out header information
        //****************************************
        _VTUFile << "<?xml version=\"1.0\"?>\n";
        _VTUFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
        _VTUFile << "<UnstructuredGrid>\n";
        _VTUFile << "<Piece NumberOfPoints=\"" << mesh.GetBulkMeshNodesNum() << "\" NumberOfCells=\"" << mesh.GetBulkMeshBulkElmtsNum() << "\">\n";
        _VTUFile << "<Points>\n";
        _VTUFile << "<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";

        //*****************************
        // print out node coordinates
        _VTUFile<<scientific << setprecision(6);
        for (i = 1; i <= mesh.GetBulkMeshNodesNum(); ++i){
            _VTUFile << mesh.GetBulkMeshIthNodeJthCoord(i, 1) << " ";
            _VTUFile << mesh.GetBulkMeshIthNodeJthCoord(i, 2) << " ";
            _VTUFile << mesh.GetBulkMeshIthNodeJthCoord(i, 3) << "\n";
        }
        _VTUFile << "</DataArray>\n";
        _VTUFile << "</Points>\n";

        //***************************************
        //*** For cell information
        //***************************************
        _VTUFile << "<Cells>\n";
        _VTUFile << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        for (e = 1; e <= mesh.GetBulkMeshBulkElmtsNum(); ++e){
            for (j = 1; j <= mesh.GetBulkMeshIthBulkElmtNodesNum(e); ++j){
                _VTUFile << mesh.GetBulkMeshIthBulkElmtJthNodeID(e, j) - 1 << " ";
            }
            _VTUFile << "\n";
        }
        _VTUFile << "</DataArray>\n";

        // for offset
        _VTUFile << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        int offset = 0;
        for (e = 1; e <= mesh.GetBulkMeshBulkElmtsNum(); ++e){
            offset += mesh.GetBulkMeshIthBulkElmtNodesNum(e);
            _VTUFile << offset << "\n";
        }
        _VTUFile << "</DataArray>\n";

        // For connectivity
        _VTUFile << "<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for (e = 1; e <= mesh.GetBulkMeshBulkElmtsNum(); ++e){
            _VTUFile << mesh.GetBulkMeshIthBulkElmtVTKCellType(e) << "\n";
        }
        _VTUFile << "</DataArray>\n";
        _VTUFile << "</Cells>\n";

        // for point data
        _VTUFile << "<PointData Scalar=\"sol\">\n";

        string dofname;
        PetscScalar value;

        for (j = 1; j <= dofHandler.GetDofsNumPerNode(); ++j){
            dofname = dofHandler.GetIthDofName(j);
            _VTUFile << "<DataArray type=\"Float64\" Name=\"" << dofname << "\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            _VTUFile << scientific << setprecision(6);
            for (i = 1; i <= mesh.GetBulkMeshNodesNum(); ++i){
                iInd = dofHandler.GetIthNodeJthDofIndex(i, j) - 1;
                VecGetValues(_Useq,1,&iInd,&value);
                _VTUFile << value << "\n";
            }
            _VTUFile << "</DataArray>\n\n";
        }

        //************************************
        //*** for projection value
        //************************************
        for(j=1;j<=nProj;++j){
            _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<projname[j-1]<<"\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            for(i=1;i<=mesh.GetBulkMeshNodesNum();++i){
                iInd =(i-1)*(1+nProj)+j;
                VecGetValues(_ProjSeq,1,&iInd,&value);
                _VTUFile<<scientific<<setprecision(6)<<value<<"\n";
            }
            _VTUFile<<"</DataArray>\n\n";
        }

        _VTUFile << "</PointData>\n";
        _VTUFile << "</Piece>\n";
        _VTUFile << "</UnstructuredGrid>\n";
        _VTUFile << "</VTKFile>" << endl;

        _VTUFile.close();
    }

    VecScatterDestroy(&_scatterU);
    VecDestroy(&_Useq);

    VecScatterDestroy(&_scatterProj);
    VecDestroy(&_ProjSeq);
}