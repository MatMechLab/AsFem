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

void OutputSystem::WriteResult2VTU(const Mesh &mesh,const DofHandler &dofHandler,const SolutionSystem &solutionSystem){
    MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);

    VecScatterCreateToAll(solutionSystem._Unew,&_scatterU,&_Useq);
    VecScatterBegin(_scatterU,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterU,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);


    //*** for projected variables
    VecScatterCreateToAll(solutionSystem._Proj,&_scatterProj,&_ProjSeq);
    VecScatterBegin(_scatterProj,solutionSystem._Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterProj,solutionSystem._Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);

    //*** for projected scalar materials
    VecScatterCreateToAll(solutionSystem._ProjScalarMate,&_scatterProjScalar,&_ProjScalarSeq);
    VecScatterBegin(_scatterProjScalar,solutionSystem._ProjScalarMate,_ProjScalarSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterProjScalar,solutionSystem._ProjScalarMate,_ProjScalarSeq,INSERT_VALUES,SCATTER_FORWARD);

    //*** for projected vector materials
    VecScatterCreateToAll(solutionSystem._ProjVectorMate,&_scatterProjVector,&_ProjVectorSeq);
    VecScatterBegin(_scatterProjVector,solutionSystem._ProjVectorMate,_ProjVectorSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterProjVector,solutionSystem._ProjVectorMate,_ProjVectorSeq,INSERT_VALUES,SCATTER_FORWARD);

    //*** for projected rank-2 materials
    VecScatterCreateToAll(solutionSystem._ProjRank2Mate,&_scatterProjRank2,&_ProjRank2Seq);
    VecScatterBegin(_scatterProjRank2,solutionSystem._ProjRank2Mate,_ProjRank2Seq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterProjRank2,solutionSystem._ProjRank2Mate,_ProjRank2Seq,INSERT_VALUES,SCATTER_FORWARD);

    //*** for projected rank-4 materials
    VecScatterCreateToAll(solutionSystem._ProjRank4Mate,&_scatterProjRank4,&_ProjRank4Seq);
    VecScatterBegin(_scatterProjRank4,solutionSystem._ProjRank4Mate,_ProjRank4Seq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterProjRank4,solutionSystem._ProjRank4Mate,_ProjRank4Seq,INSERT_VALUES,SCATTER_FORWARD);

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
        _VTUFile << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n";
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

        // for our solutions and projected quantities
        string ScalarName,VectorName,Rank2Name,Rank4Name,TensorName;

        ScalarName="<PointData Scalar=\"";
        for(i=1;i<=dofHandler.GetDofsNumPerNode();i++){
            ScalarName+=dofHandler.GetIthDofName(i)+" ";
        }
        // for projected variables
        if(solutionSystem.GetProjNumPerNode()>0){
            for(auto it:solutionSystem.GetProjNameVec()){
                ScalarName+=it+" ";
            }
        }
        // for projected scalar materials
        if(solutionSystem.GetScalarMateProjNumPerNode()>0){
            for(auto it:solutionSystem.GetScalarMateNameVec()){
                ScalarName+=it+" ";
            }
        }
        ScalarName+="\" ";

        // now we do the same procedure for vector materials
        if(solutionSystem.GetVectorMateProjNumPerNode()<1){
            VectorName="";
        }
        else{
            VectorName="Vector=\"";
            for(auto it:solutionSystem.GetVectorMateNameVec()){
                VectorName+=it+" ";
            }
            VectorName+="\" ";
        }
        TensorName="";
        if(solutionSystem.GetRank2MateProjNumPerNode()<1){
            Rank2Name="";
        }
        else{
            TensorName="Tensor=";
            Rank2Name="";
            for(auto it:solutionSystem.GetRank2MateNameVec()){
                Rank2Name+=it+" ";
            }
        }
        if(solutionSystem.GetRank4MateProjNumPerNode()<1){
            Rank4Name="";
        }
        else{
            TensorName="Tensor=";
            Rank4Name="";
            for(auto it:solutionSystem.GetRank4MateNameVec()){
                Rank4Name+=it+" ";
            }
        }
        if(solutionSystem.GetRank2MateProjNumPerNode()||solutionSystem.GetRank4MateProjNumPerNode()){
            TensorName="Tensor=\"";
            TensorName+=Rank2Name+Rank4Name+"\" ";
        }

        _VTUFile<<ScalarName<<VectorName<<TensorName<<">\n";

        string dofname,projname;
        PetscScalar value;
        int nProj;

        // output solutions
        for (j = 1;j<=dofHandler.GetDofsNumPerNode();++j){
            dofname =dofHandler.GetIthDofName(j);
            _VTUFile<<"<DataArray type=\"Float64\" Name=\"" << dofname << "\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            _VTUFile<<scientific<<setprecision(6);
            for (i = 1; i <= mesh.GetBulkMeshNodesNum(); ++i){
                iInd = dofHandler.GetBulkMeshIthNodeJthDofIndex(i, j) - 1;
                VecGetValues(_Useq,1,&iInd,&value);
                _VTUFile << value << "\n";
            }
            _VTUFile << "</DataArray>\n\n";
        }

        //************************************
        //*** for projected variables
        //************************************
        nProj=solutionSystem.GetProjNumPerNode();
        for(j=1;j<=nProj;++j){
            projname=solutionSystem.GetIthProjName(j);
            _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<projname<<"\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            for(i=1;i<=mesh.GetBulkMeshNodesNum();++i){
                iInd =(i-1)*(1+nProj)+j;
                VecGetValues(_ProjSeq,1,&iInd,&value);
                _VTUFile<<scientific<<setprecision(6)<<value<<"\n";
            }
            _VTUFile<<"</DataArray>\n\n";
        }
        //************************************
        //*** for projected scalar variables
        //************************************
        nProj=solutionSystem.GetScalarMateProjNumPerNode();
        for(j=1;j<=nProj;++j){
            projname=solutionSystem.GetIthScalarMateName(j);
            _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<projname<<"\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            for(i=1;i<=mesh.GetBulkMeshNodesNum();++i){
                iInd =(i-1)*(1+nProj)+j;
                VecGetValues(_ProjScalarSeq,1,&iInd,&value);
                _VTUFile<<scientific<<setprecision(6)<<value<<"\n";
            }
            _VTUFile<<"</DataArray>\n\n";
        }
        //************************************
        //*** for projected vector variables
        //************************************
        nProj=solutionSystem.GetVectorMateProjNumPerNode();
        for(j=1;j<=nProj;++j){
            projname=solutionSystem.GetIthVectorMateName(j);
            _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<projname<<"\"  NumberOfComponents=\"3\" format=\"ascii\">\n";
            for(i=1;i<=mesh.GetBulkMeshNodesNum();++i){
                iInd =(i-1)*(1+nProj*3)+3*(j-1)+1;
                VecGetValues(_ProjVectorSeq,1,&iInd,&value);
                _VTUFile<<scientific<<setprecision(6)<<value<<" ";
                //***************************************
                iInd =(i-1)*(1+nProj*3)+3*(j-1)+2;
                VecGetValues(_ProjVectorSeq,1,&iInd,&value);
                _VTUFile<<scientific<<setprecision(6)<<value<<" ";
                //***************************************
                iInd =(i-1)*(1+nProj*3)+3*(j-1)+3;
                VecGetValues(_ProjVectorSeq,1,&iInd,&value);
                _VTUFile<<scientific<<setprecision(6)<<value<<"\n";
            }
            _VTUFile<<"</DataArray>\n\n";
        }
        //************************************
        //*** for projected rank-2 tensor variables
        //************************************
        nProj=solutionSystem.GetRank2MateProjNumPerNode();
        int i1,j1,ii;
        for(j=1;j<=nProj;++j){
            projname=solutionSystem.GetIthRank2MateName(j);
            _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<projname<<"\"  NumberOfComponents=\"9\" format=\"ascii\">\n";
            for(i=1;i<=mesh.GetBulkMeshNodesNum();++i){
                ii=0;
                for(i1=1;i1<=3;i1++){
                    for(j1=1;j1<=3;j1++){
                        ii+=1;
                        iInd =(i-1)*(1+nProj*9)+9*(j-1)+ii;
                        VecGetValues(_ProjRank2Seq,1,&iInd,&value);
                        _VTUFile<<scientific<<setprecision(6)<<value<<" ";
                    }
                    _VTUFile<<"\n";
                }
                _VTUFile<<"\n";
            }
            _VTUFile<<"</DataArray>\n\n";
        }
        //************************************
        //*** for projected rank-4 tensor variables
        //************************************
        nProj=solutionSystem.GetRank4MateProjNumPerNode();
        for(j=1;j<=nProj;++j){
            projname=solutionSystem.GetIthRank4MateName(j);
            _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<projname<<"\"  NumberOfComponents=\"36\" format=\"ascii\">\n";
            for(i=1;i<=mesh.GetBulkMeshNodesNum();++i){
                ii=0;
                for(i1=1;i1<=6;i1++){
                    for(j1=1;j1<=6;j1++){
                        ii+=1;
                        iInd =(i-1)*(1+nProj*36)+36*(j-1)+ii;
                        VecGetValues(_ProjRank4Seq,1,&iInd,&value);
                        _VTUFile<<scientific<<setprecision(6)<<value<<" ";
                    }
                    _VTUFile<<"\n";
                }
                _VTUFile<<"\n";
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
    // for projected variables
    VecScatterDestroy(&_scatterProj);
    VecDestroy(&_ProjSeq);
    // for projected scalars
    VecScatterDestroy(&_scatterProjScalar);
    VecDestroy(&_ProjScalarSeq);
    // for projected vectors
    VecScatterDestroy(&_scatterProjVector);
    VecDestroy(&_ProjVectorSeq);
    // for projected rank-2
    VecScatterDestroy(&_scatterProjRank2);
    VecDestroy(&_ProjRank2Seq);
    // for projected rank-4
    VecScatterDestroy(&_scatterProjRank4);
    VecDestroy(&_ProjRank4Seq);
}
void OutputSystem::WriteResult2VTU(const int &step, const Mesh &mesh, const DofHandler &dofHandler,const SolutionSystem &solutionSystem){
    MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);

    VecScatterCreateToAll(solutionSystem._Unew,&_scatterU,&_Useq);
    VecScatterBegin(_scatterU,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterU,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);

    //*** for projected variables
    VecScatterCreateToAll(solutionSystem._Proj,&_scatterProj,&_ProjSeq);
    VecScatterBegin(_scatterProj,solutionSystem._Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterProj,solutionSystem._Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);

    //*** for projected scalar materials
    VecScatterCreateToAll(solutionSystem._ProjScalarMate,&_scatterProjScalar,&_ProjScalarSeq);
    VecScatterBegin(_scatterProjScalar,solutionSystem._ProjScalarMate,_ProjScalarSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterProjScalar,solutionSystem._ProjScalarMate,_ProjScalarSeq,INSERT_VALUES,SCATTER_FORWARD);

    //*** for projected vector materials
    VecScatterCreateToAll(solutionSystem._ProjVectorMate,&_scatterProjVector,&_ProjVectorSeq);
    VecScatterBegin(_scatterProjVector,solutionSystem._ProjVectorMate,_ProjVectorSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterProjVector,solutionSystem._ProjVectorMate,_ProjVectorSeq,INSERT_VALUES,SCATTER_FORWARD);

    //*** for projected rank-2 materials
    VecScatterCreateToAll(solutionSystem._ProjRank2Mate,&_scatterProjRank2,&_ProjRank2Seq);
    VecScatterBegin(_scatterProjRank2,solutionSystem._ProjRank2Mate,_ProjRank2Seq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterProjRank2,solutionSystem._ProjRank2Mate,_ProjRank2Seq,INSERT_VALUES,SCATTER_FORWARD);

    //*** for projected rank-4 materials
    VecScatterCreateToAll(solutionSystem._ProjRank4Mate,&_scatterProjRank4,&_ProjRank4Seq);
    VecScatterBegin(_scatterProjRank4,solutionSystem._ProjRank4Mate,_ProjRank4Seq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterProjRank4,solutionSystem._ProjRank4Mate,_ProjRank4Seq,INSERT_VALUES,SCATTER_FORWARD);

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
        _VTUFile << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\">\n";
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

        // for our solutions and projected quantities
        string ScalarName,VectorName,Rank2Name,Rank4Name,TensorName;

        ScalarName="<PointData Scalar=\"";
        for(i=1;i<=dofHandler.GetDofsNumPerNode();i++){
            ScalarName+=dofHandler.GetIthDofName(i)+" ";
        }
        // for projected variables
        if(solutionSystem.GetProjNumPerNode()>0){
            for(auto it:solutionSystem.GetProjNameVec()){
                ScalarName+=it+" ";
            }
        }
        // for projected scalar materials
        if(solutionSystem.GetScalarMateProjNumPerNode()>0){
            for(auto it:solutionSystem.GetScalarMateNameVec()){
                ScalarName+=it+" ";
            }
        }
        ScalarName+="\" ";

        // now we do the same procedure for vector materials
        if(solutionSystem.GetVectorMateProjNumPerNode()<1){
            VectorName="";
        }
        else{
            VectorName="Vector=\"";
            for(auto it:solutionSystem.GetVectorMateNameVec()){
                VectorName+=it+" ";
            }
            VectorName+="\" ";
        }
        TensorName="";
        if(solutionSystem.GetRank2MateProjNumPerNode()<1){
            Rank2Name="";
        }
        else{
            TensorName="Tensor=";
            Rank2Name="";
            for(auto it:solutionSystem.GetRank2MateNameVec()){
                Rank2Name+=it+" ";
            }
        }
        if(solutionSystem.GetRank4MateProjNumPerNode()<1){
            Rank4Name="";
        }
        else{
            TensorName="Tensor=";
            Rank4Name="";
            for(auto it:solutionSystem.GetRank4MateNameVec()){
                Rank4Name+=it+" ";
            }
        }
        if(solutionSystem.GetRank2MateProjNumPerNode()||solutionSystem.GetRank4MateProjNumPerNode()){
            TensorName="Tensor=\"";
            TensorName+=Rank2Name+Rank4Name+"\" ";
        }

        _VTUFile<<ScalarName<<VectorName<<TensorName<<">\n";

        string dofname,projname;
        PetscScalar value;
        int nProj;

        // output solutions
        for (j = 1;j<=dofHandler.GetDofsNumPerNode();++j){
            dofname =dofHandler.GetIthDofName(j);
            _VTUFile<<"<DataArray type=\"Float64\" Name=\"" << dofname << "\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            _VTUFile<<scientific<<setprecision(6);
            for (i = 1; i <= mesh.GetBulkMeshNodesNum(); ++i){
                iInd = dofHandler.GetBulkMeshIthNodeJthDofIndex(i, j) - 1;
                VecGetValues(_Useq,1,&iInd,&value);
                _VTUFile << value << "\n";
            }
            _VTUFile << "</DataArray>\n\n";
        }

        //************************************
        //*** for projected variables
        //************************************
        nProj=solutionSystem.GetProjNumPerNode();
        for(j=1;j<=nProj;++j){
            projname=solutionSystem.GetIthProjName(j);
            _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<projname<<"\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            for(i=1;i<=mesh.GetBulkMeshNodesNum();++i){
                iInd =(i-1)*(1+nProj)+j;
                VecGetValues(_ProjSeq,1,&iInd,&value);
                _VTUFile<<scientific<<setprecision(6)<<value<<"\n";
            }
            _VTUFile<<"</DataArray>\n\n";
        }
        //************************************
        //*** for projected scalar variables
        //************************************
        nProj=solutionSystem.GetScalarMateProjNumPerNode();
        for(j=1;j<=nProj;++j){
            projname=solutionSystem.GetIthScalarMateName(j);
            _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<projname<<"\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
            for(i=1;i<=mesh.GetBulkMeshNodesNum();++i){
                iInd =(i-1)*(1+nProj)+j;
                VecGetValues(_ProjScalarSeq,1,&iInd,&value);
                _VTUFile<<scientific<<setprecision(6)<<value<<"\n";
            }
            _VTUFile<<"</DataArray>\n\n";
        }
        //************************************
        //*** for projected vector variables
        //************************************
        nProj=solutionSystem.GetVectorMateProjNumPerNode();
        for(j=1;j<=nProj;++j){
            projname=solutionSystem.GetIthVectorMateName(j);
            _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<projname<<"\"  NumberOfComponents=\"3\" format=\"ascii\">\n";
            for(i=1;i<=mesh.GetBulkMeshNodesNum();++i){
                iInd =(i-1)*(1+nProj*3)+3*(j-1)+1;
                VecGetValues(_ProjVectorSeq,1,&iInd,&value);
                _VTUFile<<scientific<<setprecision(6)<<value<<" ";
                //***************************************
                iInd =(i-1)*(1+nProj*3)+3*(j-1)+2;
                VecGetValues(_ProjVectorSeq,1,&iInd,&value);
                _VTUFile<<scientific<<setprecision(6)<<value<<" ";
                //***************************************
                iInd =(i-1)*(1+nProj*3)+3*(j-1)+3;
                VecGetValues(_ProjVectorSeq,1,&iInd,&value);
                _VTUFile<<scientific<<setprecision(6)<<value<<"\n";
            }
            _VTUFile<<"</DataArray>\n\n";
        }
        //************************************
        //*** for projected rank-2 tensor variables
        //************************************
        nProj=solutionSystem.GetRank2MateProjNumPerNode();
        int i1,j1,ii;
        for(j=1;j<=nProj;++j){
            projname=solutionSystem.GetIthRank2MateName(j);
            _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<projname<<"\"  NumberOfComponents=\"9\" format=\"ascii\">\n";
            for(i=1;i<=mesh.GetBulkMeshNodesNum();++i){
                ii=0;
                for(i1=1;i1<=3;i1++){
                    for(j1=1;j1<=3;j1++){
                        ii+=1;
                        iInd =(i-1)*(1+nProj*9)+9*(j-1)+ii;
                        VecGetValues(_ProjRank2Seq,1,&iInd,&value);
                        _VTUFile<<scientific<<setprecision(6)<<value<<" ";
                    }
                    _VTUFile<<"\n";
                }
                _VTUFile<<"\n";
            }
            _VTUFile<<"</DataArray>\n\n";
        }
        //************************************
        //*** for projected rank-4 tensor variables
        //************************************
        nProj=solutionSystem.GetRank4MateProjNumPerNode();
        for(j=1;j<=nProj;++j){
            projname=solutionSystem.GetIthRank4MateName(j);
            _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<projname<<"\"  NumberOfComponents=\"36\" format=\"ascii\">\n";
            for(i=1;i<=mesh.GetBulkMeshNodesNum();++i){
                ii=0;
                for(i1=1;i1<=6;i1++){
                    for(j1=1;j1<=6;j1++){
                        ii+=1;
                        iInd =(i-1)*(1+nProj*36)+36*(j-1)+ii;
                        VecGetValues(_ProjRank4Seq,1,&iInd,&value);
                        _VTUFile<<scientific<<setprecision(6)<<value<<" ";
                    }
                    _VTUFile<<"\n";
                }
                _VTUFile<<"\n";
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
    // for projected variables
    VecScatterDestroy(&_scatterProj);
    VecDestroy(&_ProjSeq);
    // for projected scalars
    VecScatterDestroy(&_scatterProjScalar);
    VecDestroy(&_ProjScalarSeq);
    // for projected vectors
    VecScatterDestroy(&_scatterProjVector);
    VecDestroy(&_ProjVectorSeq);
    // for projected rank-2
    VecScatterDestroy(&_scatterProjRank2);
    VecDestroy(&_ProjRank2Seq);
    // for projected rank-4
    VecScatterDestroy(&_scatterProjRank4);
    VecDestroy(&_ProjRank4Seq);
}
