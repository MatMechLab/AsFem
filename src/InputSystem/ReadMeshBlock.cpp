//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.30
//+++ Purpose: This function can read the [mesh] block from our
//+++          input file.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadMeshBlock(ifstream &in,string str,int &linenum,Mesh &mesh){
    // mesh block format:
    // 1. built-in mesh
    // [mesh]
    //   type=asfem // this line must be the first one!!!
    //   dim=2
    //   nx=10
    //   ny=10
    //   xmin=0.0
    //   xmax=1.0
    //   ymin=0.0
    //   ymax=1.0
    //   meshtype=quad4
    //   printmesh=true[default value is false]
    //   savemesh=false[default ]
    //  [end]
    // 2. use gmsh file
    //  [mesh]
    //    type=gmsh
    //    file=gmsh.msh
    //  [end]
    //********************************************************************************
    double xmin=0.0,xmax=1.0;
    double ymin=0.0,ymax=1.0;
    double zmin=0.0,zmax=1.0;
    int dim=2;
    int nx=2,ny=2,nz=2;
    bool IsSuccess=false;
    vector<double> numbers;
    bool HasXmin=false,HasXmax=false;
    bool HasYmin=false,HasYmax=false;
    bool HasZmin=false,HasZmax=false;
    bool HasDim=false;
    bool HasNx=false,HasNy=false,HasNz=false;
    bool IsBuiltIn=false;
    bool IsSaveMesh=false;
    string meshtypename;
    bool HasType=false;
    bool IsPrint=false,IsDepPrint=false;
    char buff[55];

    // now str should contain "[mesh]"
    getline(in,str);linenum+=1;
    str=StringUtils::RemoveStrSpace(str);
    str=StringUtils::StrToLower(str);
    if(str.find("type=asfem")!=string::npos){
        HasType=true;
        IsBuiltIn=true;
        getline(in,str);linenum+=1;
        str=StringUtils::StrToLower(str);
        // read the format for built-in
        while(str.find("[end]")==string::npos){
            str=StringUtils::RemoveStrSpace(str);
            str=StringUtils::StrToLower(str);
            if(StringUtils::IsCommentLine(str) || str.length()<1) {
                getline(in,str);linenum+=1;
                continue;
            }
            if(str.find("dim=")!=string::npos){
                numbers=StringUtils::SplitStrNum(str);
                if(numbers.size()<1){
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("dim=1 [2,3] should be given in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                else{
                    dim=int(numbers[0]);
                    if(dim<1||dim>3){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("invalid dim value, dim=1[2,3] should be given in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
                    else{
                        HasDim=true;
                    } 
                }
            }
            else if(str.find("xmin=")!=string::npos||
                    str.find("XMIN=")!=string::npos||
                    str.find("Xmin=")!=string::npos){
                numbers=StringUtils::SplitStrNum(str);
                if(numbers.size()<1){
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("xmin= real value should be given in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                else{
                    xmin=numbers[0];
                    HasXmin=true;
                }
            }
            else if(str.find("xmax=")!=string::npos||
                    str.find("XMAX=")!=string::npos||
                    str.find("Xmax=")!=string::npos){
                numbers=StringUtils::SplitStrNum(str);
                if(numbers.size()<1){
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("xmax= real value should be given in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                else{
                    xmax=numbers[0];
                    HasXmax=true;
                }
            }
            else if(str.find("ymin=")!=string::npos||
                    str.find("YMIN=")!=string::npos||
                    str.find("Ymin=")!=string::npos){
                numbers=StringUtils::SplitStrNum(str);
                if(numbers.size()<1){
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("ymin= real value should be given in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                else{
                    ymin=numbers[0];
                    HasYmin=true;
                }
            }
            else if(str.find("ymax=")!=string::npos||
                    str.find("YMAX=")!=string::npos||
                    str.find("Ymax=")!=string::npos){
                numbers=StringUtils::SplitStrNum(str);
                if(numbers.size()<1){
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("ymax= real value should be given in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                else{
                    ymax=numbers[0];
                    HasYmax=true;
                }
            }
            else if(str.find("zmin=")!=string::npos||
                    str.find("ZMIN=")!=string::npos||
                    str.find("Zmin=")!=string::npos){
                numbers=StringUtils::SplitStrNum(str);
                if(numbers.size()<1){
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("zmin= real value should be given in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                else{
                    zmin=numbers[0];
                    HasZmin=true;
                }
            }
            else if(str.find("zmax=")!=string::npos||
                    str.find("ZMAX=")!=string::npos||
                    str.find("Zmax=")!=string::npos){
                numbers=StringUtils::SplitStrNum(str);
                if(numbers.size()<1){
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("zmax= real value should be given in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                else{
                    zmax=numbers[0];
                    HasZmax=true;
                }
            }
            else if(str.find("nx=")!=string::npos||
                    str.find("NX=")!=string::npos||
                    str.find("Nx=")!=string::npos){
                numbers=StringUtils::SplitStrNum(str);
                if(numbers.size()<1){
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("nx= integer value should be given in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                else{
                    nx=int(numbers[0]);
                    if(nx<1){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("nx= is invalid in [mesh] block, nx= integer value should be given in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
                    HasNx=true;
                }
            }
            else if(str.find("ny=")!=string::npos||
                    str.find("NY=")!=string::npos||
                    str.find("Ny=")!=string::npos){
                numbers=StringUtils::SplitStrNum(str);
                if(numbers.size()<1){
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("ny= integer value should be given in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                else{
                    ny=int(numbers[0]);
                    if(ny<1){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("ny= is invalid in [mesh] block, ny= integer value should be given in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
                    HasNy=true;
                }
            }
            else if(str.find("nz=")!=string::npos||
                    str.find("NZ=")!=string::npos||
                    str.find("Nz=")!=string::npos){
                numbers=StringUtils::SplitStrNum(str);
                if(numbers.size()<1){
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("nz= integer value should be given in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                else{
                    nz=int(numbers[0]);
                    if(nz<1){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("nz= is invalid in [mesh] block, nz= integer value should be given in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
                    HasNz=true;
                }
            }
            else if(str.find("meshtype=")!=string::npos||
                    str.find("MESHTYPE=")!=string::npos||
                    str.find("MeshType=")!=string::npos||
                    str.find("Meshtype=")!=string::npos){
                if(str.length()<13){
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("invalid meshtyp= information in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                else{
                    meshtypename=str.substr(9,str.length());
                }
            }
            else if(str.find("savemesh=")!=string::npos||
                    str.find("Savemesh=")!=string::npos||
                    str.find("SaveMesh=")!=string::npos||
                    str.find("SAVEMESH=")!=string::npos){
                int i=str.find_first_of('=');
                string substr=str.substr(i+1,str.length());
                substr=StringUtils::RemoveStrSpace(substr);
                if(substr.find("true")!=string::npos||
                   substr.find("True")!=string::npos||
                   substr.find("TRUE")!=string::npos){
                    IsSaveMesh=true;
                }
                else if(substr.find("false")!=string::npos||
                        substr.find("False")!=string::npos||
                        substr.find("FALSE")!=string::npos){
                    IsSaveMesh=false;
                }
                else{
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("unsupported option in savemesh= in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
            }
            else if(str.find("printmesh=")!=string::npos||
                    str.find("Printmesh=")!=string::npos||
                    str.find("PrintMesh=")!=string::npos||
                    str.find("PRINTMESH=")!=string::npos){
                int i=str.find_first_of('=');
                string substr=str.substr(i+1,str.length());
                substr=StringUtils::RemoveStrSpace(substr);
                if(substr.find("true")!=string::npos||
                   substr.find("True")!=string::npos||
                   substr.find("TRUE")!=string::npos){
                    IsPrint=true;IsDepPrint=false;
                }
                else if(substr.find("false")!=string::npos||
                        substr.find("False")!=string::npos||
                        substr.find("FALSE")!=string::npos){
                    IsPrint=false;IsDepPrint=false;
                }
                else if(substr.find("dep")!=string::npos||
                        substr.find("Dep")!=string::npos||
                        substr.find("DEP")!=string::npos){
                    IsPrint=true;IsDepPrint=true;
                }
                else{
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("unsupported option in printmesh= in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
            }
            else if(str.find("[end]")!=string::npos||str.find("[END]")!=string::npos){
                break;
            }
            else{
                snprintf(buff,55,"line-%d has some errors",linenum);
                MessagePrinter::PrintErrorTxt(string(buff));
                MessagePrinter::PrintErrorTxt("unknown option in the [mesh] block");
                MessagePrinter::AsFem_Exit();
            }
            getline(in,str);linenum+=1;
        }

        //*****************************************************
        if(!HasDim){
            MessagePrinter::PrintErrorTxt("dim=1[2,3] should be given in the [mesh] block");
            MessagePrinter::AsFem_Exit();
        }

        if(dim==1){
            if(!HasNx){
                MessagePrinter::PrintErrorTxt("nx should be given for dim=1 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasXmin&&xmin!=0.0){
                MessagePrinter::PrintErrorTxt("xmin=val should be given for dim=1 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasXmax&&xmax!=1.0){
                MessagePrinter::PrintErrorTxt("xmax=val should be given for dim=1 case");
                MessagePrinter::AsFem_Exit();
            }
            if(HasNy||HasNz){
                MessagePrinter::PrintErrorTxt("for dim=1 case, you only need nx");
                MessagePrinter::AsFem_Exit();
            }
            if(HasYmin||HasYmax||HasZmin||HasZmax){
                MessagePrinter::PrintErrorTxt("for dim=1 case, you only need xmin,xmax");
                MessagePrinter::AsFem_Exit();
            }

            if(meshtypename!="edge2"&&meshtypename!="edge3"&&meshtypename!="edge4"){
                MessagePrinter::PrintErrorTxt("unsupported 1d mesh type in the [mesh] block, meshtype=edge2,edge3,edge4 is expected");
                MessagePrinter::AsFem_Exit();
            }
            IsSuccess=true;
            mesh.SetDim(dim);
            mesh.SetNx(nx);
            mesh.SetXmin(xmin);
            mesh.SetXmax(xmax);
            mesh.SetMeshTypeName(meshtypename);
        }
        else if(dim==2){
            if(!HasNx){
                MessagePrinter::PrintErrorTxt("nx should be given for dim=2 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasNy){
                MessagePrinter::PrintErrorTxt("ny should be given for dim=2 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasXmin&&xmin!=0.0){
                MessagePrinter::PrintErrorTxt("xmin=val should be given for dim=1 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasXmax&&xmax!=1.0){
                MessagePrinter::PrintErrorTxt("xmax=val should be given for dim=1 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasYmin&&ymin!=0.0) {
                MessagePrinter::PrintErrorTxt("ymin=val should be given for dim=1 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasYmax&&ymax!=1.0){
                MessagePrinter::PrintErrorTxt("ymax=val should be given for dim=2 case");
                MessagePrinter::AsFem_Exit();
            }
            if(HasNz){
                MessagePrinter::PrintErrorTxt("for dim=2 case, you only need nx and ny");
                MessagePrinter::AsFem_Exit();
            }
            if(HasZmin||HasZmax){
                MessagePrinter::PrintErrorTxt("for dim=2 case, you only need xmin/xmax and ymin/ymax");
                MessagePrinter::AsFem_Exit();
            }
            if(meshtypename!="quad4"&&meshtypename!="quad8"&&meshtypename!="quad9"){
                MessagePrinter::PrintErrorTxt("unsupported 2d mesh type in the [mesh] block, meshtype=quad4,quad8,quad9 is expected");
                MessagePrinter::AsFem_Exit();
            }
            IsSuccess=true;
            mesh.SetDim(dim);
            mesh.SetNx(nx);
            mesh.SetNy(ny);
            mesh.SetXmin(xmin);
            mesh.SetXmax(xmax);
            mesh.SetYmin(ymin);
            mesh.SetYmax(ymax);
            mesh.SetMeshTypeName(meshtypename);
        }
        else if(dim==3){
            if(!HasNx){
                MessagePrinter::PrintErrorTxt("nx should be given for dim=3 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasNy){
                MessagePrinter::PrintErrorTxt("nx should be given for dim=3 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasNy){
                MessagePrinter::PrintErrorTxt("nz should be given for dim=3 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasXmin&&xmin!=0.0){
                MessagePrinter::PrintErrorTxt("xmin=val should be given for dim=3 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasXmax&&xmax!=1.0){
                MessagePrinter::PrintErrorTxt("xmax=val should be given for dim=3 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasYmin&&ymin!=0.0){
                MessagePrinter::PrintErrorTxt("ymin=val should be given for dim=3 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasYmax&&ymax!=1.0){
                MessagePrinter::PrintErrorTxt("ymax=val should be given for dim=3 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasZmin&&zmin!=0.0){
                MessagePrinter::PrintErrorTxt("zmin=val should be given for dim=3 case");
                MessagePrinter::AsFem_Exit();
            }
            if(!HasZmax&&zmax!=1.0){
                MessagePrinter::PrintErrorTxt("zmax=val should be given for dim=3 case");
                MessagePrinter::AsFem_Exit();
            }

            if(meshtypename!="hex8"&&meshtypename!="hex20"&&meshtypename!="hex27"){
                MessagePrinter::PrintErrorTxt("unsupported 3d mesh type in the [mesh] block, meshtype=hex8,hex20,hex27 is expected");
                MessagePrinter::AsFem_Exit();
            }
            IsSuccess=true;
            mesh.SetDim(dim);
            mesh.SetNx(nx);
            mesh.SetNy(ny);
            mesh.SetNz(nz);
            mesh.SetXmin(xmin);
            mesh.SetXmax(xmax);
            mesh.SetYmin(ymin);
            mesh.SetYmax(ymax);
            mesh.SetZmin(zmin);
            mesh.SetZmax(zmax);
            mesh.SetMeshTypeName(meshtypename);
        }
        IsBuiltIn=true;
        IsSuccess=false;
    }
    else if(str.find("type=gmsh")!=string::npos||
            str.find("type=Gmsh")!=string::npos||
            str.find("type=GMSH")!=string::npos){
        HasType=true;
        IsBuiltIn=false;
        getline(in,str);linenum+=1;
        str=StringUtils::StrToLower(str);
        bool HasFileName=false;
        while(str.find("[end]")==string::npos&&
              str.find("[END]")==string::npos){
            str=StringUtils::RemoveStrSpace(str);
            if(StringUtils::IsCommentLine(str)||str.length()<1){
                getline(in,str);linenum+=1;
                str=StringUtils::StrToLower(str);
                continue;
            }
            if(str.find("file=")!=string::npos||
               str.find("FILE=")!=string::npos){
                if(str.compare(str.length()-4,4,".msh")==0||
                   str.compare(str.length()-4,4,".Msh")==0||
                   str.compare(str.length()-4,4,".MSH")==0){
                    string filename=str.substr(5,str.length());
                    _meshio.SetMeshFileName(filename);
                    HasFileName=true;
                    MessagePrinter::PrintNormalTxt("   start to import mesh ...");
                    IsSuccess=_meshio.ReadMeshFromFile(mesh);
                    MessagePrinter::PrintNormalTxt("   import mesh finished   !");
                    // mesh.SetMeshMode(false);
                }
            }
            else if(str.find("savemesh=")!=string::npos||
                    str.find("Savemesh=")!=string::npos||
                    str.find("SaveMesh=")!=string::npos||
                    str.find("SAVEMESH=")!=string::npos){
                int i=str.find_first_of('=');
                string substr=str.substr(i+1,str.length());
                substr=StringUtils::RemoveStrSpace(substr);
                if(substr.find("true")!=string::npos||
                   substr.find("True")!=string::npos||
                   substr.find("TRUE")!=string::npos){
                    IsSaveMesh=true;
                }
                else if(substr.find("false")!=string::npos||
                        substr.find("False")!=string::npos||
                        substr.find("FALSE")!=string::npos){
                    IsSaveMesh=false;
                }
                else{
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("unsupported option in savemesh= in the [mesh], savemesh=true[false] is expected");
                    MessagePrinter::AsFem_Exit();
                }
            }
            else if(str.find("printmesh=")!=string::npos||
                    str.find("Printmesh=")!=string::npos||
                    str.find("PrintMesh=")!=string::npos||
                    str.find("PRINTMESH=")!=string::npos){
                int i=str.find_first_of('=');
                string substr=str.substr(i+1,str.length());
                substr=StringUtils::RemoveStrSpace(substr);
                if(substr.find("true")!=string::npos||
                   substr.find("True")!=string::npos||
                   substr.find("TRUE")!=string::npos){
                    IsPrint=true;IsDepPrint=false;
                }
                else if(substr.find("false")!=string::npos||
                        substr.find("False")!=string::npos||
                        substr.find("FALSE")!=string::npos){
                    IsPrint=false;IsDepPrint=false;
                }
                else if(substr.find("dep")!=string::npos||
                        substr.find("Dep")!=string::npos||
                        substr.find("DEP")!=string::npos){
                    IsPrint=true;IsDepPrint=true;
                }
                else{
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("unsupported option in printmesh= in the [mesh] block, printmesh=true[false] is expected");
                    MessagePrinter::AsFem_Exit();
                }
            }
            else if(str.find("[end]")!=string::npos||
                    str.find("[END]")!=string::npos){
                break;
            }
            else if(str.find("[]")!=string::npos){
                snprintf(buff,55,"line-%d has some errors",linenum);
                MessagePrinter::PrintErrorTxt(string(buff));
                MessagePrinter::PrintErrorTxt("bracket pair is not complete in the [mesh] block");
                MessagePrinter::AsFem_Exit();
            }
            else{
                snprintf(buff,55,"line-%d has some errors",linenum);
                MessagePrinter::PrintErrorTxt(string(buff));
                MessagePrinter::PrintErrorTxt("unknown option in the [mesh] block");
                MessagePrinter::AsFem_Exit();
            }
            getline(in,str);linenum+=1;
        }
        if(!HasFileName){
            IsSuccess=false;
            snprintf(buff,55,"line-%d has some errors",linenum);
            MessagePrinter::PrintErrorTxt(string(buff));
            MessagePrinter::PrintErrorTxt("file=correct file name should be given in the [mesh] block");
            MessagePrinter::AsFem_Exit();
        }
    }
    else if(str.find("type=abaqus")!=string::npos||
            str.find("type=Abaqus")!=string::npos||
            str.find("type=ABAQUS")!=string::npos){
        HasType=true;
        IsBuiltIn=false;
        getline(in,str);linenum+=1;
        str=StringUtils::StrToLower(str);
        bool HasFileName=false;
        while(str.find("[end]")==string::npos&&
              str.find("[END]")==string::npos){
            str=StringUtils::RemoveStrSpace(str);
            if(StringUtils::IsCommentLine(str)||str.length()<1){
                getline(in,str);linenum+=1;
                str=StringUtils::StrToLower(str);
                continue;
            }
            if(str.find("file=")!=string::npos||
               str.find("FILE=")!=string::npos){
                if(str.compare(str.length()-4,4,".inp")==0||
                   str.compare(str.length()-4,4,".Inp")==0||
                   str.compare(str.length()-4,4,".INP")==0){
                    string filename=str.substr(5,str.length());
                    _meshio.SetMeshFileName(filename);
                    MessagePrinter::PrintNormalTxt("   start to import mesh ...");
                    IsSuccess=_meshio.ReadMeshFromFile(mesh);
                    MessagePrinter::PrintNormalTxt("   import mesh finished   !");
                    HasFileName=true;
                }
            }
            else if(str.find("savemesh=")!=string::npos||
                    str.find("Savemesh=")!=string::npos||
                    str.find("SaveMesh=")!=string::npos||
                    str.find("SAVEMESH=")!=string::npos){
                int i=str.find_first_of('=');
                string substr=str.substr(i+1,str.length());
                substr=StringUtils::RemoveStrSpace(substr);
                if(substr.find("true")!=string::npos||
                   substr.find("True")!=string::npos||
                   substr.find("TRUE")!=string::npos){
                    IsSaveMesh=true;
                }
                else if(substr.find("false")!=string::npos||
                        substr.find("False")!=string::npos||
                        substr.find("FALSE")!=string::npos){
                    IsSaveMesh=false;
                }
                else{
                    IsSuccess=false;
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("unsupported option in savemesh= in the [mesh] block, savemesh=true[false] is expected");
                    MessagePrinter::AsFem_Exit();
                }
            }
            else if(str.find("printmesh=")!=string::npos||
                    str.find("Printmesh=")!=string::npos||
                    str.find("PrintMesh=")!=string::npos||
                    str.find("PRINTMESH=")!=string::npos){
                int i=str.find_first_of('=');
                string substr=str.substr(i+1,str.length());
                substr=StringUtils::RemoveStrSpace(substr);
                if(substr.find("true")!=string::npos||
                   substr.find("True")!=string::npos||
                   substr.find("TRUE")!=string::npos){
                    IsPrint=true;IsDepPrint=false;
                }
                else if(substr.find("false")!=string::npos||
                        substr.find("False")!=string::npos||
                        substr.find("FALSE")!=string::npos){
                    IsPrint=false;IsDepPrint=false;
                }
                else if(substr.find("dep")!=string::npos||
                        substr.find("Dep")!=string::npos||
                        substr.find("DEP")!=string::npos){
                    IsPrint=true;IsDepPrint=true;
                }
                else{
                    snprintf(buff,55,"line-%d has some errors",linenum);
                    MessagePrinter::PrintErrorTxt(string(buff));
                    MessagePrinter::PrintErrorTxt("unsupported option in printmesh= in the [mesh] block, printmesh=true[false] is expected");
                    MessagePrinter::AsFem_Exit();
                }
            }
            else if(str.find("[end]")!=string::npos||
                    str.find("[END]")!=string::npos){
                break;
            }
            else if(str.find("[]")!=string::npos){
                snprintf(buff,55,"line-%d has some errors",linenum);
                MessagePrinter::PrintErrorTxt(string(buff));
                MessagePrinter::PrintErrorTxt("block bracket pair is not complete in the [mesh] block");
                MessagePrinter::AsFem_Exit();
            }
            else{
                snprintf(buff,55,"line-%d has some errors",linenum);
                MessagePrinter::PrintErrorTxt(string(buff));
                MessagePrinter::PrintErrorTxt("unknown option in the [mesh] block");
                MessagePrinter::AsFem_Exit();
            }
            getline(in,str);linenum+=1;
        }
        if(!HasFileName){
            IsSuccess=false;
            MessagePrinter::PrintErrorTxt("file=correct file name should be given in the [mesh] block");
            MessagePrinter::AsFem_Exit();
        }
    }


    string substr=_InputFileName.substr(0,_InputFileName.find_first_of('.'))+"_mesh.vtu";
    _MeshFileName=substr;

    if(!HasType){
        MessagePrinter::PrintErrorTxt("no type= found in the [mesh] block, type=asfem[gmsh] should be given");
        MessagePrinter::AsFem_Exit();
    }

    
    if(IsBuiltIn){
        MessagePrinter::PrintNormalTxt("   start to create mesh ...");
        if(!mesh.CreateMesh()){
            MessagePrinter::PrintErrorTxt("create mesh failed !!! please check your input file");
            MessagePrinter::AsFem_Exit();
        }
        IsSuccess=true;
        MessagePrinter::PrintNormalTxt("   mesh generation finished !");
    }
    if(IsSaveMesh){
        mesh.SaveLagrangeMesh(_MeshFileName);
        snprintf(buff,55,"save mesh to [%39s]",substr.c_str());
        MessagePrinter::PrintNormalTxt(string(buff));   
    }
    if(IsPrint){
        if(IsDepPrint){
            mesh.PrintMeshDetailInfo();
        }
        else{
            mesh.PrintMeshInfo();
        }
    }

    return IsSuccess;
}