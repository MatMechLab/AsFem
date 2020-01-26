#include "OutputSystem/OutputSystem.h"

void OutputSystem::WriteResultToVTU(Mesh &mesh,DofHandler &dofHandler,const Eigen::VectorXd &U){
    _VTUFileName=_OutputFilePrefix+".vtu";
    _VTUFile.open(_VTUFileName,ios::out);
    if(!_VTUFile.is_open()){
        cout<<"*** Error: can\'t create a new vtu file(="<<_VTUFileName<<")!!!"<<endl;
        cout<<"***        please make sure you have write permission!***"<<endl;
        Msg_AsFem_Exit();
    }
    int i,j,iInd,e;

    //****************************************
    //*** print out header information
    //****************************************
    _VTUFile<<"<?xml version=\"1.0\"?>\n";
    _VTUFile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
    _VTUFile<<"<UnstructuredGrid>\n";
    _VTUFile<<"<Piece NumberOfPoints=\""<<mesh.GetNodesNum()<<"\" NumberOfCells=\""<<mesh.GetBulkElmtsNum()<<"\">\n";
    _VTUFile<<"<Points>\n";
    _VTUFile<<"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";

    //*****************************
    // print out node coordinates
    _VTUFile<<scientific<<setprecision(6);
    for(i=1;i<=mesh.GetNodesNum();++i){
        _VTUFile<<mesh.GetIthNodeJthCoord(i,1)<<" ";
        _VTUFile<<mesh.GetIthNodeJthCoord(i,2)<<" ";
        _VTUFile<<mesh.GetIthNodeJthCoord(i,3)<<"\n";
    }
    _VTUFile<<"</DataArray>\n";
    _VTUFile<<"</Points>\n";


    //***************************************
    //*** For cell information
    //***************************************
    _VTUFile<<"<Cells>\n";
    _VTUFile<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        for(j=1;j<=mesh.GetIthBulkElmtNodesNum(e);++j){
            _VTUFile<<mesh.GetIthBulkElmtJthConn(e,j)-1<<" ";
        }
        _VTUFile<<"\n";
    }
    _VTUFile<<"</DataArray>\n";

    // for offset
    _VTUFile<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    int offset=0;
    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        offset+=mesh.GetIthBulkElmtNodesNum(e);
        _VTUFile<<offset<<"\n";
    }
    _VTUFile<<"</DataArray>\n";

    // For connectivity
    _VTUFile<<"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        _VTUFile<<mesh.GetIthBulkElmtVTKType(e)<<"\n";
    }
     _VTUFile<<"</DataArray>\n";
    _VTUFile<<"</Cells>\n";

    // for point data
    _VTUFile<<"<PointData Scalar=\"sol\">\n";

    string dofname;
    
    for(j=1;j<=dofHandler.GetDofsNumPerNode();++j){
        dofname=dofHandler.GetIthDofNameFromList(j);
        _VTUFile<<"<DataArray type=\"Float64\" Name=\""<<dofname<<"\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
        _VTUFile<<scientific<<setprecision(6);
        for(i=1;i<=mesh.GetNodesNum();++i){
            iInd=dofHandler.GetIthNodeJthDofIndex(i,j)-1;
            _VTUFile<<U(iInd)<<"\n";
        }
        _VTUFile<<"</DataArray>\n\n";
    }
    _VTUFile<<"</PointData>\n";
    _VTUFile<<"</Piece>\n";
    _VTUFile<<"</UnstructuredGrid>\n";
    _VTUFile<<"</VTKFile>"<<endl;

    _VTUFile.close();

}
//*** for output with projection(static case)
void OutputSystem::WriteResultToVTU(Mesh &mesh,DofHandler &dofHandler,const Eigen::VectorXd &U,const int &nproj,const vector<string> &namelist,const Eigen::MatrixXd &Proj){
    _VTUFileName=_OutputFilePrefix+".vtu";
    _VTUFile.open(_VTUFileName,ios::out);
    if(!_VTUFile.is_open()){
        cout<<"*** Error: can\'t create a new vtu file(="<<_VTUFileName<<")!!!"<<endl;
        cout<<"***        please make sure you have write permission!***"<<endl;
        Msg_AsFem_Exit();
    }
    int i,j,iInd,e;

    //****************************************
    //*** print out header information
    //****************************************
    _VTUFile<<"<?xml version=\"1.0\"?>\n";
    _VTUFile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
    _VTUFile<<"<UnstructuredGrid>\n";
    _VTUFile<<"<Piece NumberOfPoints=\""<<mesh.GetNodesNum()<<"\" NumberOfCells=\""<<mesh.GetBulkElmtsNum()<<"\">\n";
    _VTUFile<<"<Points>\n";
    _VTUFile<<"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";;

    //*****************************
    // print out node coordinates
    _VTUFile<<scientific<<setprecision(6);
    for(i=1;i<=mesh.GetNodesNum();++i){
        _VTUFile<<mesh.GetIthNodeJthCoord(i,1)<<" ";
        _VTUFile<<mesh.GetIthNodeJthCoord(i,2)<<" ";
        _VTUFile<<mesh.GetIthNodeJthCoord(i,3)<<"\n";
    }
    _VTUFile<<"</DataArray>\n";
    _VTUFile<<"</Points>\n";


    //***************************************
    //*** For cell information
    //***************************************
    _VTUFile<<"<Cells>\n";
    _VTUFile<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        for(j=1;j<=mesh.GetIthBulkElmtNodesNum(e);++j){
            _VTUFile<<mesh.GetIthBulkElmtJthConn(e,j)-1<<" ";
        }
        _VTUFile<<"\n";
    }
    _VTUFile<<"</DataArray>\n";

    // for offset
    _VTUFile<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    int offset=0;
    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        offset+=mesh.GetIthBulkElmtNodesNum(e);
        _VTUFile<<offset<<"\n";
    }
    _VTUFile<<"</DataArray>\n";

    // For connectivity
    _VTUFile<<"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        _VTUFile<<mesh.GetIthBulkElmtVTKType(e)<<"\n";
    }
     _VTUFile<<"</DataArray>\n";
    _VTUFile<<"</Cells>\n";

    // for point data
    _VTUFile<<"<PointData Scalar=\"sol\">\n";

    string dofname;
    
    for(j=1;j<=dofHandler.GetDofsNumPerNode();++j){
        dofname=dofHandler.GetIthDofNameFromList(j);
        _VTUFile<<"<DataArray type=\"Float64\" Name=\""<<dofname<<"\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
        _VTUFile<<scientific<<setprecision(6);
        for(i=1;i<=mesh.GetNodesNum();++i){
            iInd=dofHandler.GetIthNodeJthDofIndex(i,j)-1;
            _VTUFile<<U(iInd)<<"\n";
        }
        _VTUFile<<"</DataArray>\n\n";
    }
    //************************************
    //*** for projection value
    //************************************
    double value,weight;
    for(j=1;j<=nproj;++j){
        _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<namelist[j-1]<<"\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
        for(i=1;i<=mesh.GetNodesNum();++i){
            value=Proj.coeff(i-1,j);
            weight=Proj.coeff(i-1,0);
            // if(weight){};
            // this method give the value smaller than real one!
            // if(abs(value/weight)>1.0e-13){
            //     _VTUFile<<scientific<<setprecision(6)<<value/weight<<"\n";
            // }
            // value*weight gives even smaller values than correct one!!!
            // if(abs(value*weight)>1.0e-13){
            //     _VTUFile<<scientific<<setprecision(6)<<value*weight<<"\n";
            // }
            // this value is not correct
            // if(abs(value)>1.0e-13){
            //     _VTUFile<<scientific<<setprecision(6)<<value<<"\n";
            // }
            // value/weight is the closest one to the correct value, anyhow
            // it still has problems!!!
            if(abs(value/weight)>1.0e-13){
                _VTUFile<<scientific<<setprecision(6)<<value/weight<<"\n";
            }
            else{
                _VTUFile<<scientific<<setprecision(6)<<1.0e-13<<"\n";
            }
        }
        _VTUFile<<"</DataArray>\n\n";
    }


    // for the end part
    _VTUFile<<"</PointData>\n";
    _VTUFile<<"</Piece>\n";
    _VTUFile<<"</UnstructuredGrid>\n";
    _VTUFile<<"</VTKFile>"<<endl;

    _VTUFile.close();
}
//*******************************************
void OutputSystem::WriteResultToVTU(const int &step,Mesh &mesh,DofHandler &dofHandler,const Eigen::VectorXd &U){
    
    ostringstream ss;
    ss<<setfill('0')<<setw(8)<<step;
    _VTUFileName=_OutputFilePrefix+"-"+ss.str()+".vtu";
    _VTUFile.open(_VTUFileName,ios::out);
    if(!_VTUFile.is_open()){
        cout<<"*** Error: can\'t create a new vtu file(="<<_VTUFileName<<")!!!"<<endl;
        cout<<"***        please make sure you have write permission!***"<<endl;
        Msg_AsFem_Exit();
    }
    int i,j,iInd,e;

    //****************************************
    //*** print out header information
    //****************************************
    _VTUFile<<"<?xml version=\"1.0\"?>\n";
    _VTUFile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
    _VTUFile<<"<UnstructuredGrid>\n";
    _VTUFile<<"<Piece NumberOfPoints=\""<<mesh.GetNodesNum()<<"\" NumberOfCells=\""<<mesh.GetBulkElmtsNum()<<"\">\n";
    _VTUFile<<"<Points>\n";
    _VTUFile<<"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";

    //*****************************
    // print out node coordinates
    _VTUFile<<scientific<<setprecision(6);
    for(i=1;i<=mesh.GetNodesNum();++i){
        _VTUFile<<mesh.GetIthNodeJthCoord(i,1)<<" ";
        _VTUFile<<mesh.GetIthNodeJthCoord(i,2)<<" ";
        _VTUFile<<mesh.GetIthNodeJthCoord(i,3)<<"\n";
    }
    _VTUFile<<"</DataArray>\n";
    _VTUFile<<"</Points>\n";


    //***************************************
    //*** For cell information
    //***************************************
    _VTUFile<<"<Cells>\n";
    _VTUFile<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        for(j=1;j<=mesh.GetIthBulkElmtNodesNum(e);++j){
            _VTUFile<<mesh.GetIthBulkElmtJthConn(e,j)-1<<" ";
        }
        _VTUFile<<"\n";
    }
    _VTUFile<<"</DataArray>\n";

    // for offset
    _VTUFile<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    int offset=0;
    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        offset+=mesh.GetIthBulkElmtNodesNum(e);
        _VTUFile<<offset<<"\n";
    }
    _VTUFile<<"</DataArray>\n";

    // For connectivity
    _VTUFile<<"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        _VTUFile<<mesh.GetIthBulkElmtVTKType(e)<<"\n";
    }
     _VTUFile<<"</DataArray>\n";
    _VTUFile<<"</Cells>\n";

    // for point data
    _VTUFile<<"<PointData Scalar=\"sol\">\n";

    string dofname;
    
    for(j=1;j<=dofHandler.GetDofsNumPerNode();++j){
        dofname=dofHandler.GetIthDofNameFromList(j);
        _VTUFile<<"<DataArray type=\"Float64\" Name=\""<<dofname<<"\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
        _VTUFile<<scientific<<setprecision(6);
        for(i=1;i<=mesh.GetNodesNum();++i){
            iInd=dofHandler.GetIthNodeJthDofIndex(i,j)-1;
            _VTUFile<<U(iInd)<<"\n";
        }
        _VTUFile<<"</DataArray>\n\n";
    }
    _VTUFile<<"</PointData>\n";
    _VTUFile<<"</Piece>\n";
    _VTUFile<<"</UnstructuredGrid>\n";
    _VTUFile<<"</VTKFile>"<<endl;

    _VTUFile.close();

}
//**** for output with projection in transient case
//*** for output with projection(static case)
void OutputSystem::WriteResultToVTU(const int &step,Mesh &mesh,DofHandler &dofHandler,const Eigen::VectorXd &U,const int &nproj,const vector<string> &namelist,const Eigen::MatrixXd &Proj){
    ostringstream ss;
    ss<<setfill('0')<<setw(8)<<step;
    _VTUFileName=_OutputFilePrefix+"-"+ss.str()+".vtu";
    _VTUFile.open(_VTUFileName,ios::out);
    if(!_VTUFile.is_open()){
        cout<<"*** Error: can\'t create a new vtu file(="<<_VTUFileName<<")!!!"<<endl;
        cout<<"***        please make sure you have write permission!***"<<endl;
        Msg_AsFem_Exit();
    }
    int i,j,iInd,e;

    //****************************************
    //*** print out header information
    //****************************************
    _VTUFile<<"<?xml version=\"1.0\"?>\n";
    _VTUFile<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
    _VTUFile<<"<UnstructuredGrid>\n";
    _VTUFile<<"<Piece NumberOfPoints=\""<<mesh.GetNodesNum()<<"\" NumberOfCells=\""<<mesh.GetBulkElmtsNum()<<"\">\n";
    _VTUFile<<"<Points>\n";
    _VTUFile<<"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";

    //*****************************
    // print out node coordinates
    _VTUFile<<scientific<<setprecision(6);
    for(i=1;i<=mesh.GetNodesNum();++i){
        _VTUFile<<mesh.GetIthNodeJthCoord(i,1)<<" ";
        _VTUFile<<mesh.GetIthNodeJthCoord(i,2)<<" ";
        _VTUFile<<mesh.GetIthNodeJthCoord(i,3)<<"\n";
    }
    _VTUFile<<"</DataArray>\n";
    _VTUFile<<"</Points>\n";


    //***************************************
    //*** For cell information
    //***************************************
    _VTUFile<<"<Cells>\n";
    _VTUFile<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        for(j=1;j<=mesh.GetIthBulkElmtNodesNum(e);++j){
            _VTUFile<<mesh.GetIthBulkElmtJthConn(e,j)-1<<" ";
        }
        _VTUFile<<"\n";
    }
    _VTUFile<<"</DataArray>\n";

    // for offset
    _VTUFile<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    int offset=0;
    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        offset+=mesh.GetIthBulkElmtNodesNum(e);
        _VTUFile<<offset<<"\n";
    }
    _VTUFile<<"</DataArray>\n";

    // For connectivity
    _VTUFile<<"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        _VTUFile<<mesh.GetIthBulkElmtVTKType(e)<<"\n";
    }
     _VTUFile<<"</DataArray>\n";
    _VTUFile<<"</Cells>\n";

    // for point data
    _VTUFile<<"<PointData Scalar=\"sol\">\n";

    string dofname;
    
    for(j=1;j<=dofHandler.GetDofsNumPerNode();++j){
        dofname=dofHandler.GetIthDofNameFromList(j);
        _VTUFile<<"<DataArray type=\"Float64\" Name=\""<<dofname<<"\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
        _VTUFile<<scientific<<setprecision(6);
        for(i=1;i<=mesh.GetNodesNum();++i){
            iInd=dofHandler.GetIthNodeJthDofIndex(i,j)-1;
            _VTUFile<<U(iInd)<<"\n";
        }
        _VTUFile<<"</DataArray>\n\n";
    }
    //************************************
    //*** for projection value
    //************************************
    double value,weight;
    for(j=1;j<=nproj;++j){
        _VTUFile<<"<DataArray type=\"Float64\"  Name=\""<<namelist[j-1]<<"\"  NumberOfComponents=\"1\" format=\"ascii\">\n";
        for(i=1;i<=mesh.GetNodesNum();++i){
            value=Proj.coeff(i-1,j);
            weight=Proj.coeff(i-1,0);
            // if(weight){};
            // this method give the value smaller than real one!
            // if(abs(value/weight)>1.0e-13){
            //     _VTUFile<<scientific<<setprecision(6)<<value/weight<<"\n";
            // }
            // value*weight gives even smaller values than correct one!!!
            // if(abs(value*weight)>1.0e-13){
            //     _VTUFile<<scientific<<setprecision(6)<<value*weight<<"\n";
            // }
            // this value is not correct
            // if(abs(value)>1.0e-13){
            //     _VTUFile<<scientific<<setprecision(6)<<value<<"\n";
            // }
            // value/weight is the closest one to the correct value, anyhow
            // it still has problems!!!
            if(abs(value/weight)>1.0e-13){
                _VTUFile<<scientific<<setprecision(6)<<value/weight<<"\n";
            }
            else{
                _VTUFile<<scientific<<setprecision(6)<<1.0e-13<<"\n";
            }
        }
        _VTUFile<<"</DataArray>\n\n";
    }


    // for the end part
    _VTUFile<<"</PointData>\n";
    _VTUFile<<"</Piece>\n";
    _VTUFile<<"</UnstructuredGrid>\n";
    _VTUFile<<"</VTKFile>"<<endl;

    _VTUFile.close();
}