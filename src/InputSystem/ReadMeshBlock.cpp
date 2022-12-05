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
//+++ Date   : 2022.05.07
//+++ Purpose: the input file reading system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"
#include "Mesh/MeshGenerator.h"
#include "Mesh/MeshFileImporter.h"

bool InputSystem::readMeshBlock(nlohmann::json &t_json,Mesh &t_mesh){
    // the json already read 'mesh' !!!
    MeshType meshtype;
    int dim,nx,ny,nz;
    double xmin,xmax,ymin,ymax,zmin,zmax;
    bool HasNx,HasNy,HasNz,HasMeshType,IsSaveMesh;
    string meshfile;

    nx=2;ny=2;nz=2;
    xmin=ymin=zmin=0.0;
    xmax=ymax=zmax=1.0;
    dim=-1;
    meshtype=MeshType::NULLTYPE;
    HasNx=false;HasNy=false;HasNz=false;
    HasMeshType=false;IsSaveMesh=false;
    if(t_json.contains("type")){
        if(!t_json.at("type").is_string()){
            MessagePrinter::printErrorTxt("the type name of your mesh block is not a valid string");
            return false;
        }
        string meshtypename=t_json.at("type");
        if(meshtypename.find("asfem")!=string::npos){
            // for AsFem's built-in type
            if(t_json.contains("dim")){
                if(!t_json.at("dim").is_number_integer()){
                    MessagePrinter::printErrorTxt("the dim number in your mesh block is not a valid integer");
                    return false;
                }
                dim=static_cast<int>(t_json.at("dim"));
                if(dim==1){
                    // 1d case
                    if(t_json.contains("nx")){
                        if(!t_json.at("nx").is_number_integer()){
                            MessagePrinter::printErrorTxt("the nx number in your mesh block is not a valid integer");
                            return false;
                        }
                        nx=static_cast<int>(t_json.at("nx"));
                        HasNx=true;
                    }
                    if(t_json.contains("xmin")){
                        if(!t_json.at("xmin").is_number_float()){
                            MessagePrinter::printErrorTxt("the xmin value in your mesh block is not a valid float");
                            return false;
                        }
                        xmin=static_cast<double>(t_json.at("xmin"));
                    }
                    else{
                        xmin=0.0;
                    }

                    if(t_json.contains("xmax")){
                        if(!t_json.at("xmax").is_number_float()){
                            MessagePrinter::printErrorTxt("the xmax value in your mesh block is not a valid float");
                            return false;
                        }
                        xmax=static_cast<double>(t_json.at("xmax"));
                    }
                    else{
                        xmax=1.0;
                    }
                    if(t_json.contains("meshtype")){
                        if(!t_json.at("meshtype").is_string()){
                            MessagePrinter::printErrorTxt("the meshtype in your mesh block is not a valid string");
                            return false;
                        }
                        string meshname=t_json.at("meshtype");
                        if(meshname=="edge2"){
                            meshtype=MeshType::EDGE2;
                            HasMeshType=true;
                        }
                        else if(meshname=="edge3"){
                            meshtype=MeshType::EDGE3;
                            HasMeshType=true;
                        }
                        else if(meshname=="edge4"){
                            meshtype=MeshType::EDGE4;
                            HasMeshType=true;
                        }
                        else if(meshname=="edge5"){
                            meshtype=MeshType::EDGE5;
                            HasMeshType=true;
                        }
                        else{
                            MessagePrinter::printErrorTxt("unsupported meshtype in your mesh block, it should be edge2,edge3,edge4");
                            return false;
                        }
                    }// end-of-meshtype-read

                    if(t_json.contains("savemesh")){
                        if(!t_json.at("savemesh").is_boolean()){
                            MessagePrinter::printErrorTxt("invalid boolean value for savemesh in your mesh block, it should be true/false");
                            return false;
                        }
                        IsSaveMesh=static_cast<bool>(t_json.at("savemesh"));
                    }
                    // now all the necessary information is ready for 1d case, we start to check unnecessary inputs
                    if(t_json.contains("ny")){
                        MessagePrinter::printErrorTxt("ny is invalid for 1d case in your mesh block, please check your input file");
                        return false;
                    }
                    if(t_json.contains("nz")){
                        MessagePrinter::printErrorTxt("nz is invalid for 1d case in your mesh block, please check your input file");
                        return false;
                    }

                    if(t_json.contains("ymin")){
                        MessagePrinter::printErrorTxt("ymin is invalid for 1d case in your mesh block, please check your input file");
                        return false;
                    }
                    if(t_json.contains("ymax")){
                        MessagePrinter::printErrorTxt("ymax is invalid for 1d case in your mesh block, please check your input file");
                        return false;
                    }
                     if(t_json.contains("zmin")){
                        MessagePrinter::printErrorTxt("zmin is invalid for 1d case in your mesh block, please check your input file");
                        return false;
                    }
                    if(t_json.contains("zmax")){
                        MessagePrinter::printErrorTxt("zmax is invalid for 1d case in your mesh block, please check your input file");
                        return false;
                    }
                    if(!HasNx){
                        MessagePrinter::printErrorTxt("for 1d case, nx must be given in your mesh block, please check your input file");
                        return false;
                    }
                    if(!HasMeshType){
                        MessagePrinter::printErrorTxt("for 1d case, meshtype must be given in your mesh block, please check your input file");
                        return false;
                    }
                    // assign the mesh info to mesh data
                    t_mesh.getBulkMeshMeshDataRef().m_nx=nx;
                    t_mesh.getBulkMeshMeshDataRef().m_xmin=xmin;
                    t_mesh.getBulkMeshMeshDataRef().m_xmax=xmax;
                    t_mesh.getBulkMeshMeshDataRef().m_ymin=0.0;
                    t_mesh.getBulkMeshMeshDataRef().m_ymax=0.0;
                    t_mesh.getBulkMeshMeshDataRef().m_zmin=0.0;
                    t_mesh.getBulkMeshMeshDataRef().m_zmax=0.0;
                }// end-of-dim1-case
                else if(dim==2){
                    // for 2d case
                    if(t_json.contains("nx")){
                        if(!t_json.at("nx").is_number_integer()){
                            MessagePrinter::printErrorTxt("the nx number in your mesh block is not a valid integer");
                            return false;
                        }
                        nx=static_cast<int>(t_json.at("nx"));
                        HasNx=true;
                    }
                    if(t_json.contains("ny")){
                        if(!t_json.at("ny").is_number_integer()){
                            MessagePrinter::printErrorTxt("the ny number in your mesh block is not a valid integer");
                            return false;
                        }
                        ny=static_cast<int>(t_json.at("ny"));
                        HasNy=true;
                    }
                    if(t_json.contains("xmin")){
                        if(!t_json.at("xmin").is_number_float()){
                            MessagePrinter::printErrorTxt("the xmin value in your mesh block is not a valid float");
                            return false;
                        }
                        xmin=static_cast<double>(t_json.at("xmin"));
                    }
                    else{
                        xmin=0.0;
                    }
                    if(t_json.contains("xmax")){
                        if(!t_json.at("xmax").is_number_float()){
                            MessagePrinter::printErrorTxt("the xmax value in your mesh block is not a valid float");
                            return false;
                        }
                        xmax=static_cast<double>(t_json.at("xmax"));
                    }
                    else{
                        xmax=1.0;
                    }
                    if(t_json.contains("ymin")){
                        if(!t_json.at("ymin").is_number_float()){
                            MessagePrinter::printErrorTxt("the ymin value in your mesh block is not a valid float");
                            return false;
                        }
                        ymin=static_cast<double>(t_json.at("ymin"));
                    }
                    else{
                        ymin=0.0;
                    }
                    if(t_json.contains("ymax")){
                        if(!t_json.at("ymax").is_number_float()){
                            MessagePrinter::printErrorTxt("the ymax value in your mesh block is not a valid float");
                            return false;
                        }
                        ymax=static_cast<double>(t_json.at("ymax"));
                    }
                    else{
                        ymax=1.0;
                    }
                    if(t_json.contains("meshtype")){
                        if(!t_json.at("meshtype").is_string()){
                            MessagePrinter::printErrorTxt("the meshtype in your mesh block is not a valid string");
                            return false;
                        }
                        string meshname=t_json.at("meshtype");
                        if(meshname=="quad4"){
                            meshtype=MeshType::QUAD4;
                            HasMeshType=true;
                        }
                        else if(meshname=="quad8"){
                            meshtype=MeshType::QUAD8;
                            HasMeshType=true;
                        }
                        else if(meshname=="quad9"){
                            meshtype=MeshType::QUAD9;
                            HasMeshType=true;
                        }
                        else{
                            
                            MessagePrinter::printErrorTxt("unsupported 2d meshtype in your mesh block, it should be quad4,quad8,quad9");
                            return false;
                        }
                    }// end-of-meshtype-read

                    if(t_json.contains("savemesh")){
                        if(!t_json.at("savemesh").is_boolean()){
                            MessagePrinter::printErrorTxt("invalid boolean value for savemesh in your mesh block, it should be true/false");
                            return false;
                        }
                        IsSaveMesh=static_cast<bool>(t_json.at("savemesh"));
                    }
                    // now all the necessary information is ready for 2d case, we start to check unnecessary inputs
                    if(t_json.contains("nz")){
                        MessagePrinter::printErrorTxt("nz is invalid for 2d case in your mesh block, please check your input file");
                        return false;
                    }

                    if(t_json.contains("zmin")){
                        MessagePrinter::printErrorTxt("zmin is invalid for 2d case in your mesh block, please check your input file");
                        return false;
                    }
                    if(t_json.contains("zmax")){
                        MessagePrinter::printErrorTxt("zmax is invalid for 2d case in your mesh block, please check your input file");
                        return false;
                    }
                    if(!HasNx){
                        MessagePrinter::printErrorTxt("for 2d case, nx must be given in your mesh block, please check your input file");
                        return false;
                    }
                    if(!HasNy){
                        MessagePrinter::printErrorTxt("for 2d case, ny must be given in your mesh block, please check your input file");
                        return false;
                    }
                    if(!HasMeshType){
                        MessagePrinter::printErrorTxt("for 2d case, meshtype must be given in your mesh block, please check your input file");
                        return false;
                    }
                    // assign the mesh info to mesh data
                    t_mesh.getBulkMeshMeshDataRef().m_nx=nx;
                    t_mesh.getBulkMeshMeshDataRef().m_ny=ny;
                    t_mesh.getBulkMeshMeshDataRef().m_xmin=xmin;
                    t_mesh.getBulkMeshMeshDataRef().m_xmax=xmax;
                    t_mesh.getBulkMeshMeshDataRef().m_ymin=ymin;
                    t_mesh.getBulkMeshMeshDataRef().m_ymax=ymax;
                    t_mesh.getBulkMeshMeshDataRef().m_zmin=0.0;
                    t_mesh.getBulkMeshMeshDataRef().m_zmax=0.0;

                }//end-of-dim2-case
                else if(dim==3){
                    // for 2d case
                    if(t_json.contains("nx")){
                        if(!t_json.at("nx").is_number_integer()){
                            MessagePrinter::printErrorTxt("the nx number in your mesh block is not a valid integer");
                            return false;
                        }
                        nx=static_cast<int>(t_json.at("nx"));
                        HasNx=true;
                    }
                    if(t_json.contains("ny")){
                        if(!t_json.at("ny").is_number_integer()){
                            MessagePrinter::printErrorTxt("the ny number in your mesh block is not a valid integer");
                            return false;
                        }
                        ny=static_cast<int>(t_json.at("ny"));
                        HasNy=true;
                    }
                    if(t_json.contains("nz")){
                        if(!t_json.at("nz").is_number_integer()){
                            MessagePrinter::printErrorTxt("the nz number in your mesh block is not a valid integer");
                            return false;
                        }
                        nz=static_cast<int>(t_json.at("nz"));
                        HasNz=true;
                    }
                    if(t_json.contains("xmin")){
                        if(!t_json.at("xmin").is_number_float()){
                            MessagePrinter::printErrorTxt("the xmin value in your mesh block is not a valid float");
                            return false;
                        }
                        xmin=static_cast<double>(t_json.at("xmin"));
                    }
                    else{
                        xmin=0.0;
                    }
                    if(t_json.contains("xmax")){
                        if(!t_json.at("xmax").is_number_float()){
                            MessagePrinter::printErrorTxt("the xmax value in your mesh block is not a valid float");
                            return false;
                        }
                        xmax=static_cast<double>(t_json.at("xmax"));
                    }
                    else{
                        xmax=1.0;
                    }
                    if(t_json.contains("ymin")){
                        if(!t_json.at("ymin").is_number_float()){
                            MessagePrinter::printErrorTxt("the ymin value in your mesh block is not a valid float");
                            return false;
                        }
                        ymin=static_cast<double>(t_json.at("ymin"));
                    }
                    else{
                        ymin=0.0;
                    }
                    if(t_json.contains("ymax")){
                        if(!t_json.at("ymax").is_number_float()){
                            MessagePrinter::printErrorTxt("the ymax value in your mesh block is not a valid float");
                            return false;
                        }
                        ymax=static_cast<double>(t_json.at("ymax"));
                    }
                    else{
                        ymax=1.0;
                    }
                    if(t_json.contains("zmin")){
                        if(!t_json.at("zmin").is_number_float()){
                            MessagePrinter::printErrorTxt("the zmin value in your mesh block is not a valid float");
                            return false;
                        }
                        zmin=static_cast<double>(t_json.at("zmin"));
                    }
                    else{
                        zmin=0.0;
                    }
                    if(t_json.contains("zmax")){
                        if(!t_json.at("zmax").is_number_float()){
                            MessagePrinter::printErrorTxt("the zmax value in your mesh block is not a valid float");
                            return false;
                        }
                        zmax=static_cast<double>(t_json.at("zmax"));
                    }
                    else{
                        zmax=1.0;
                    }
                    if(t_json.contains("meshtype")){
                        if(!t_json.at("meshtype").is_string()){
                            MessagePrinter::printErrorTxt("the meshtype in your mesh block is not a valid string");
                            return false;
                        }
                        string meshname=t_json.at("meshtype");
                        if(meshname=="hex8"){
                            meshtype=MeshType::HEX8;
                            HasMeshType=true;
                        }
                        else if(meshname=="hex20"){
                            meshtype=MeshType::HEX20;
                            HasMeshType=true;
                        }
                        else if(meshname=="hex27"){
                            meshtype=MeshType::HEX27;
                            HasMeshType=true;
                        }
                        else{
                            
                            MessagePrinter::printErrorTxt("unsupported 3d meshtype in your mesh block, it should be hex8,hex20,hex27");
                            return false;
                        }
                    }// end-of-meshtype-read

                    if(t_json.contains("savemesh")){
                        if(!t_json.at("savemesh").is_boolean()){
                            MessagePrinter::printErrorTxt("invalid boolean value for savemesh in your mesh block, it should be true/false");
                            return false;
                        }
                        IsSaveMesh=static_cast<bool>(t_json.at("savemesh"));
                    }
                    // now all the necessary information is ready for 3d case, we start to check unnecessary inputs
                    if(!HasNx){
                        MessagePrinter::printErrorTxt("for 3d case, nx must be given in your mesh block, please check your input file");
                        return false;
                    }
                    if(!HasNy){
                        MessagePrinter::printErrorTxt("for 3d case, ny must be given in your mesh block, please check your input file");
                        return false;
                    }
                    if(!HasNz){
                        MessagePrinter::printErrorTxt("for 3d case, nz must be given in your mesh block, please check your input file");
                        return false;
                    }
                    if(!HasMeshType){
                        MessagePrinter::printErrorTxt("for 3d case, meshtype must be given in your mesh block, please check your input file");
                        return false;
                    }
                    // assign the mesh info to mesh data
                    t_mesh.getBulkMeshMeshDataRef().m_nx=nx;
                    t_mesh.getBulkMeshMeshDataRef().m_ny=ny;
                    t_mesh.getBulkMeshMeshDataRef().m_nz=nz;
                    t_mesh.getBulkMeshMeshDataRef().m_xmin=xmin;
                    t_mesh.getBulkMeshMeshDataRef().m_xmax=xmax;
                    t_mesh.getBulkMeshMeshDataRef().m_ymin=ymin;
                    t_mesh.getBulkMeshMeshDataRef().m_ymax=ymax;
                    t_mesh.getBulkMeshMeshDataRef().m_zmin=zmin;
                    t_mesh.getBulkMeshMeshDataRef().m_zmax=zmax;
                }//end-of-dim3-case
                else{
                    MessagePrinter::printErrorTxt("dim>3 or dim<1 is invalid for mesh generation, please check your input file");
                    return false;
                }

            }// end-of-contain-dim
            MeshGenerator meshGenerator;
            if(meshGenerator.createMesh(dim,meshtype,t_mesh.getBulkMeshMeshDataRef())){
                MessagePrinter::printNormalTxt("mesh generator is done, your mesh is generated");
                if(IsSaveMesh){
                    t_mesh.saveBulkMesh2VTU(m_inputfile_name);
                    MessagePrinter::printNormalTxt("save mesh to "+m_inputfile_name.substr(0,m_inputfile_name.size()-5)+"-mesh.vtu");
                    MessagePrinter::printStars();
                }
                return true;
            }
            else{
                MessagePrinter::printErrorTxt("something is wrong with your mesh generator, please check your input file or your mesh generator");
                return false;
            }
        }// end-of-asfem-type-mesh
        else if(meshtypename=="msh2"){
            // for msh file with v2.0 format
            if(!t_json.contains("file")){
                MessagePrinter::printErrorTxt("No mesh file ('file') option is given, please check your input file");
                return false;
            }
            if(!t_json.at("file").is_string()){
                MessagePrinter::printErrorTxt("'file' option is not a valid string, please check your input file");
                return false;
            }
            meshfile=t_json.at("file");
            if(meshfile.substr(meshfile.size()-4)!=".msh"){
                MessagePrinter::printErrorTxt("'file' must be '*.msh' for msh2 type, please check your input file");
                return false;
            }

            if(t_json.contains("savemesh")){
                if(!t_json.at("savemesh").is_boolean()){
                    MessagePrinter::printErrorTxt("invalid boolean value for savemesh in your mesh block, it should be true/false");
                    return false;
                }
                IsSaveMesh=static_cast<bool>(t_json.at("savemesh"));
            }
            MeshFileImporter importer;

            if(importer.importMsh2Mesh(meshfile,t_mesh.getBulkMeshMeshDataRef())){
                if(IsSaveMesh){
                    t_mesh.saveBulkMesh2VTU(m_inputfile_name);
                    MessagePrinter::printNormalTxt("save mesh to "+m_inputfile_name.substr(0,m_inputfile_name.size()-5)+"-mesh.vtu");
                }
                return true;
            }
            else{
                return false;
            }
        }
        else if(meshtypename=="msh4"){
            // for msh file with v4.0 format
            if(!t_json.contains("file")){
                MessagePrinter::printErrorTxt("No mesh file ('file') option is given, please check your input file");
                return false;
            }
            if(!t_json.at("file").is_string()){
                MessagePrinter::printErrorTxt("'file' option is not a valid string, please check your input file");
                return false;
            }
            meshfile=t_json.at("file");
            if(meshfile.substr(meshfile.size()-4)!=".msh"){
                MessagePrinter::printErrorTxt("'file' must be '*.msh' for msh4 type, please check your input file");
                return false;
            }

            if(t_json.contains("savemesh")){
                if(!t_json.at("savemesh").is_boolean()){
                    MessagePrinter::printErrorTxt("invalid boolean value for savemesh in your mesh block, it should be true/false");
                    return false;
                }
                IsSaveMesh=static_cast<bool>(t_json.at("savemesh"));
            }
            MeshFileImporter importer;

            if(importer.importMsh4Mesh(meshfile,t_mesh.getBulkMeshMeshDataRef())){
                if(IsSaveMesh){
                    t_mesh.saveBulkMesh2VTU(m_inputfile_name);
                    MessagePrinter::printNormalTxt("save mesh to "+m_inputfile_name.substr(0,m_inputfile_name.size()-5)+"-mesh.vtu");
                }
                return true;
            }
            else{
                return false;
            }
        }
        else if(meshtypename=="gmsh2"){
            // for gmsh2 file from netgen
            if(!t_json.contains("file")){
                MessagePrinter::printErrorTxt("No mesh file ('file') option is given, please check your input file");
                return false;
            }
            if(!t_json.at("file").is_string()){
                MessagePrinter::printErrorTxt("'file' option is not a valid string, please check your input file");
                return false;
            }
            meshfile=t_json.at("file");
            if(meshfile.size()<6+1){
                MessagePrinter::printErrorTxt("file="+meshfile+" is not a valid gmsh2 file, please check your input file");
                MessagePrinter::exitAsFem();
            }
            if(meshfile.substr(meshfile.size()-6)!=".gmsh2"){
                MessagePrinter::printErrorTxt("'file' must be '*.gmsh2' for gmsh2 type, please check your input file");
                return false;
            }

            if(t_json.contains("savemesh")){
                if(!t_json.at("savemesh").is_boolean()){
                    MessagePrinter::printErrorTxt("invalid boolean value for savemesh in your mesh block, it should be true/false");
                    return false;
                }
                IsSaveMesh=static_cast<bool>(t_json.at("savemesh"));
            }
            MeshFileImporter importer;

            if(importer.importGmsh2Mesh(meshfile,t_mesh.getBulkMeshMeshDataRef())){
                if(IsSaveMesh){
                    t_mesh.saveBulkMesh2VTU(m_inputfile_name);
                    MessagePrinter::printNormalTxt("save mesh to "+m_inputfile_name.substr(0,m_inputfile_name.size()-5)+"-mesh.vtu");
                }
                return true;
            }
            else{
                return false;
            }
        }
        else{
            MessagePrinter::printErrorTxt("Unsupported mesh import type, please check your input file");
            MessagePrinter::exitAsFem();
            return false;
        }
    }
    return true;
}