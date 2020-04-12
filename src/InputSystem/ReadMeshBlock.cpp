//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadMeshBlock(ifstream &in,string str,int &linenum,Mesh &mesh)
{
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
    //  [end]
    // 2. use gmsh file
    //  [mesh]
    //    type=gmsh
    //    file=gmsh.msh
    //  [end]
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
    string meshtype;
    bool HasType=false;
    bool IsPrint=false,IsDepPrint=false;;

    // now str should contain "[mesh]"
    getline(in,str);linenum+=1;
    str=RemoveStrSpace(str);
    str=StrToLower(str);
    if(str.find("type=asfem")!=string::npos){
        HasType=true;
        IsBuiltIn=true;
        getline(in,str);linenum+=1;
        str=StrToLower(str);
        // read the format for built-in
        while(str.find("[end]")==string::npos){
            str=RemoveStrSpace(str);
            str=StrToLower(str);
            if(IsCommentLine(str) || str.length()<1) {
                getline(in,str);linenum+=1;
                continue;
            }
            if(str.find("dim=")!=string::npos){
                numbers=SplitStrNum(str);
                if(numbers.size()<1){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: dim=1 [2,3] should be given in [mesh] block !!!            ***\n");
                    Msg_AsFem_Exit();
                }
                else{
                    dim=int(numbers[0]);
                    if(dim<1||dim>3){
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid dim value, dim=1[2,3] should be given in [mesh]!!! ***\n");
                    }
                    else{
                        HasDim=true;
                    } 
                }
            }
            else if(str.find("xmin=")!=string::npos||
                    str.find("XMIN=")!=string::npos||
                    str.find("Xmin=")!=string::npos){
                numbers=SplitStrNum(str);
                if(numbers.size()<1){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: xmin= real value should be given in [mesh]                 ***\n");
                    Msg_AsFem_Exit();
                }
                else{
                    xmin=numbers[0];
                    HasXmin=true;
                }
            }
            else if(str.find("xmax=")!=string::npos||
                    str.find("XMAX=")!=string::npos||
                    str.find("Xmax=")!=string::npos){
                numbers=SplitStrNum(str);
                if(numbers.size()<1){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: xmax= real value should be given in [mesh]                 ***\n");
                    Msg_AsFem_Exit();
                }
                else{
                    xmax=numbers[0];
                    HasXmax=true;
                }
            }
            else if(str.find("ymin=")!=string::npos||
                    str.find("YMIN=")!=string::npos||
                    str.find("Ymin=")!=string::npos){
                numbers=SplitStrNum(str);
                if(numbers.size()<1){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: ymin= real value should be given in [mesh]                 ***\n");
                    Msg_AsFem_Exit();
                }
                else{
                    ymin=numbers[0];
                    HasYmin=true;
                }
            }
            else if(str.find("ymax=")!=string::npos||
                    str.find("YMAX=")!=string::npos||
                    str.find("Ymax=")!=string::npos){
                numbers=SplitStrNum(str);
                if(numbers.size()<1){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: ymax= real value should be given in [mesh]                 ***\n");
                    Msg_AsFem_Exit();
                }
                else{
                    ymax=numbers[0];
                    HasYmax=true;
                }
            }
            else if(str.find("zmin=")!=string::npos||
                    str.find("ZMIN=")!=string::npos||
                    str.find("Zmin=")!=string::npos){
                numbers=SplitStrNum(str);
                if(numbers.size()<1){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: zmin= real value should be given in [mesh]                 ***\n");
                    Msg_AsFem_Exit();
                }
                else{
                    zmin=numbers[0];
                    HasZmin=true;
                }
            }
            else if(str.find("zmax=")!=string::npos||
                    str.find("ZMAX=")!=string::npos||
                    str.find("Zmax=")!=string::npos){
                numbers=SplitStrNum(str);
                if(numbers.size()<1){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: zmax= real value should be given in [mesh]                 ***\n");
                    Msg_AsFem_Exit();
                }
                else{
                    zmax=numbers[0];
                    HasZmax=true;
                }
            }
            else if(str.find("nx=")!=string::npos||
                    str.find("NX=")!=string::npos||
                    str.find("Nx=")!=string::npos){
                numbers=SplitStrNum(str);
                if(numbers.size()<1){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: nx= integer value should be given in [mesh] !!!            ***\n");
                    Msg_AsFem_Exit();
                }
                else{
                    nx=int(numbers[0]);
                    if(nx<1){
                        Msg_Input_LineError(linenum);
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: nx= is invalid in [mesh] block              !!!            ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        nx= integer value should be given in [mesh] !!!            ***\n");
                        Msg_AsFem_Exit();
                    }
                    HasNx=true;
                }
            }
            else if(str.find("ny=")!=string::npos||
                    str.find("NY=")!=string::npos||
                    str.find("Ny=")!=string::npos){
                numbers=SplitStrNum(str);
                if(numbers.size()<1){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: ny= integer value should be given in [mesh] !!!            ***\n");Msg_AsFem_Exit();
                }
                else{
                    ny=int(numbers[0]);
                    if(ny<1){
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: ny= is invalid in [mesh] block              !!!            ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        ny= integer value should be given in [mesh] !!!            ***\n");
                        Msg_AsFem_Exit();
                    }
                    HasNy=true;
                }
            }
            else if(str.find("nz=")!=string::npos||
                    str.find("NZ=")!=string::npos||
                    str.find("Nz=")!=string::npos){
                numbers=SplitStrNum(str);
                if(numbers.size()<1){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: nz= integer value should be given in [mesh] !!!            ***\n");Msg_AsFem_Exit();
                }
                else{
                    nz=int(numbers[0]);
                    if(nz<1)
                    {
                        PetscPrintf(PETSC_COMM_WORLD,"*** Error: nz= is invalid in [mesh] block              !!!            ***\n");
                        PetscPrintf(PETSC_COMM_WORLD,"***        nz= integer value should be given in [mesh] !!!            ***\n");
                        Msg_AsFem_Exit();
                    }
                    HasNz=true;
                }
            }
            else if(str.find("meshtype=")!=string::npos||
                    str.find("MESHTYPE=")!=string::npos||
                    str.find("MeshType=")!=string::npos||
                    str.find("Meshtype=")!=string::npos){
                if(str.length()<13){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid meshtyp= information in [mesh] block !!!           ***\n");
                    Msg_AsFem_Exit();
                }
                else{
                    meshtype=str.substr(9,str.length());
                }
            }
            else if(str.find("savemesh=")!=string::npos||
                    str.find("Savemesh=")!=string::npos||
                    str.find("SaveMesh=")!=string::npos||
                    str.find("SAVEMESH=")!=string::npos){
                int i=str.find_first_of('=');
                string substr=str.substr(i+1,str.length());
                substr=RemoveStrSpace(substr);
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
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported option in savemesh= in [mesh] block      !!!   ***\n");
                    Msg_AsFem_Exit();
                }
            }
            else if(str.find("printmesh=")!=string::npos||
                    str.find("Printmesh=")!=string::npos||
                    str.find("PrintMesh=")!=string::npos||
                    str.find("PRINTMESH=")!=string::npos){
                int i=str.find_first_of('=');
                string substr=str.substr(i+1,str.length());
                substr=RemoveStrSpace(substr);
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
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported option in printmesh= in [mesh] block     !!!   ***\n");
                    Msg_AsFem_Exit();
                }
            }
            else if(str.find("[end]")!=string::npos||str.find("[END]")!=string::npos){
                break;
            }
            else{
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: unknown option in [mesh] block                       !!!   ***\n");
                Msg_AsFem_Exit();
            }
            getline(in,str);linenum+=1;
        }

        //*****************************************************
        if(!HasDim){
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: dim=1[2,3] should be given in mesh block               !!! ***\n");
            Msg_AsFem_Exit();
        }

        if(dim==1){
            if(!HasNx){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: nx should be given for dim=1 case           !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasXmin&&xmin!=0.0){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: xmin=val should be given for dim=1 case     !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasXmax&&xmax!=1.0){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: xmax=val should be given for dim=1 case     !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(HasNy||HasNz){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: for dim=1 case, you only need nx            !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(HasYmin||HasYmax||HasZmin||HasZmax){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: for dim=1 case, you only need xmin,xmax     !!!            ***\n");
                Msg_AsFem_Exit();
            }

            if(meshtype!="edge2"&&meshtype!="edge3"&&meshtype!="edge4"){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported 1d mesh type in [mesh] block    !!!            ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        meshtype=edge2,edge3,edge4 is expected      !!!            ***\n");
                Msg_AsFem_Exit();
            }
            IsSuccess=true;
            mesh.SetDim(dim);
            mesh.SetNx(nx);
            mesh.SetXmin(xmin);
            mesh.SetXmax(xmax);
            mesh.SetMeshType(meshtype);
            mesh.SetMeshMode(true);
        }
        else if(dim==2){
            if(!HasNx){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: nx should be given for dim=2 case           !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasNy){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: ny should be given for dim=2 case           !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasXmin&&xmin!=0.0){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: xmin=val should be given for dim=1 case     !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasXmax&&xmax!=1.0){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: xmax=val should be given for dim=1 case     !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasYmin&&ymin!=0.0) {
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: ymin=val should be given for dim=1 case     !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasYmax&&ymax!=1.0){
                cout<<"*** Error: ymax=val should be given for dim=2 case !!!     ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(HasNz){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: for dim=2 case, you only need nx,ny         !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(HasZmin||HasZmax){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: for dim=2 case, you only need xmin, xmax, ymin, ymax   !!! ***\n");
                Msg_AsFem_Exit();
            }
            if(meshtype!="quad4"&&meshtype!="quad8"&&meshtype!="quad9"){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported 2d mesh type in [mesh] block    !!!            ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        meshtype=quad4,quad8,quad9 is expected      !!!            ***\n");
                Msg_AsFem_Exit();
            }
            IsSuccess=true;
            mesh.SetDim(dim);
            mesh.SetNx(nx);
            mesh.SetNy(ny);
            mesh.SetXmin(xmin);
            mesh.SetXmax(xmax);
            mesh.SetYmin(ymin);
            mesh.SetYmax(ymax);
            mesh.SetMeshType(meshtype);
            mesh.SetMeshMode(true);
        }
        else if(dim==3){
            if(!HasNx){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: nx should be given for dim=3 case           !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasNy){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: ny should be given for dim=3 case           !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasNy){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: nz should be given for dim=3 case           !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasXmin&&xmin!=0.0){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: xmin=val should be given for dim=3 case     !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasXmax&&xmax!=1.0){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: xmax=val should be given for dim=3 case     !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasYmin&&ymin!=0.0){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: ymin=val should be given for dim=3 case     !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasYmax&&ymax!=1.0){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: ymax=val should be given for dim=3 case     !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasZmin&&zmin!=0.0){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: zmin=val should be given for dim=3 case     !!!            ***\n");
                Msg_AsFem_Exit();
            }
            if(!HasZmax&&zmax!=1.0){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: zmax=val should be given for dim=3 case     !!!            ***\n");
                Msg_AsFem_Exit();
            }

            if(meshtype!="hex8"&&meshtype!="hex20"&&meshtype!="hex27"){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported 3d mesh type in [mesh] block    !!!            ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        meshtype=hex8,hex20,hex27 is expected       !!!            ***\n");
                Msg_AsFem_Exit();
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
            mesh.SetMeshType(meshtype);
            mesh.SetMeshMode(true);
        }
        IsBuiltIn=true;
        mesh.SetMeshMode(IsBuiltIn);
        IsSuccess=true;
    }
    else if(str.find("type=gmsh")!=string::npos||
            str.find("type=Gmsh")!=string::npos||
            str.find("type=GMSH")!=string::npos){
        HasType=true;
        IsBuiltIn=false;
        getline(in,str);linenum+=1;
        str=StrToLower(str);
        mesh.SetMeshMode(IsBuiltIn);
        mesh.SetMeshGmshMode(true);
        bool HasFileName=false;
        while(str.find("[end]")==string::npos&&
              str.find("[END]")==string::npos){
            str=RemoveStrSpace(str);
            if(IsCommentLine(str)||str.length()<1){
                getline(in,str);linenum+=1;
                str=StrToLower(str);
                continue;
            }
            if(str.find("file=")!=string::npos||
               str.find("FILE=")!=string::npos){
                if(str.compare(str.length()-4,4,".msh")==0||
                   str.compare(str.length()-4,4,".Msh")==0||
                   str.compare(str.length()-4,4,".MSH")==0){
                    string filename=str.substr(5,str.length());
                    mesh.SetMshFileName(filename);
                    IsSuccess=true;
                    HasFileName=true;
                    mesh.SetMeshMode(false);
                }
            }
            else if(str.find("savemesh=")!=string::npos||
                    str.find("Savemesh=")!=string::npos||
                    str.find("SaveMesh=")!=string::npos||
                    str.find("SAVEMESH=")!=string::npos){
                int i=str.find_first_of('=');
                string substr=str.substr(i+1,str.length());
                substr=RemoveStrSpace(substr);
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
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported option in savemesh= in [mesh]   !!!            ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        savemesh=true[false] is expected            !!!            ***\n");
                    Msg_AsFem_Exit();
                }
            }
            else if(str.find("printmesh=")!=string::npos||
                    str.find("Printmesh=")!=string::npos||
                    str.find("PrintMesh=")!=string::npos||
                    str.find("PRINTMESH=")!=string::npos){
                int i=str.find_first_of('=');
                string substr=str.substr(i+1,str.length());
                substr=RemoveStrSpace(substr);
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
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported option in printmesh= in [mesh] block     !!!   ***\n");
                    Msg_AsFem_Exit();
                }
            }
            else if(str.find("[end]")!=string::npos||
                    str.find("[END]")!=string::npos){
                break;
            }
            else if(str.find("[]")!=string::npos){
                Msg_Input_LineError(linenum);
                Msg_Input_BlockBracketNotComplete();
                Msg_AsFem_Exit();
            }
            else{
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: unknown option in [mesh] block                         !!! ***\n");
                Msg_AsFem_Exit();
            }
            getline(in,str);linenum+=1;
        }
        if(!HasFileName){
            IsSuccess=false;
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: file=correct file name should be given in [mesh] block !!! ***\n");
            Msg_AsFem_Exit();
        }
    }
    else if(str.find("type=abaqus")!=string::npos||
            str.find("type=Abaqus")!=string::npos||
            str.find("type=ABAQUS")!=string::npos){
        HasType=true;
        IsBuiltIn=false;
        getline(in,str);linenum+=1;
        str=StrToLower(str);
        mesh.SetMeshMode(IsBuiltIn);
        mesh.SetMeshAbaqusMode(true);
        bool HasFileName=false;
        while(str.find("[end]")==string::npos&&
              str.find("[END]")==string::npos){
            str=RemoveStrSpace(str);
            if(IsCommentLine(str)||str.length()<1){
                getline(in,str);linenum+=1;
                str=StrToLower(str);
                continue;
            }
            if(str.find("file=")!=string::npos||
               str.find("FILE=")!=string::npos){
                if(str.compare(str.length()-4,4,".inp")==0||
                   str.compare(str.length()-4,4,".Inp")==0||
                   str.compare(str.length()-4,4,".INP")==0){
                    string filename=str.substr(5,str.length());
                    mesh.SetMshFileName(filename);
                    IsSuccess=true;
                    HasFileName=true;
                    mesh.SetMeshMode(false);
                }
            }
            else if(str.find("savemesh=")!=string::npos||
                    str.find("Savemesh=")!=string::npos||
                    str.find("SaveMesh=")!=string::npos||
                    str.find("SAVEMESH=")!=string::npos){
                int i=str.find_first_of('=');
                string substr=str.substr(i+1,str.length());
                substr=RemoveStrSpace(substr);
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
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported option in savemesh= in [mesh]   !!!            ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        savemesh=true[false] is expected            !!!            ***\n");
                    Msg_AsFem_Exit();
                }
            }
            else if(str.find("printmesh=")!=string::npos||
                    str.find("Printmesh=")!=string::npos||
                    str.find("PrintMesh=")!=string::npos||
                    str.find("PRINTMESH=")!=string::npos){
                int i=str.find_first_of('=');
                string substr=str.substr(i+1,str.length());
                substr=RemoveStrSpace(substr);
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
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported option in printmesh= in [mesh] block     !!!   ***\n");
                    Msg_AsFem_Exit();
                }
            }
            else if(str.find("[end]")!=string::npos||
                    str.find("[END]")!=string::npos){
                break;
            }
            else if(str.find("[]")!=string::npos){
                Msg_Input_LineError(linenum);
                Msg_Input_BlockBracketNotComplete();
                Msg_AsFem_Exit();
            }
            else{
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: unknown option in [mesh] block                         !!! ***\n");
                Msg_AsFem_Exit();
            }
            getline(in,str);linenum+=1;
        }
        if(!HasFileName){
            IsSuccess=false;
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: file=correct file name should be given in [mesh] block !!! ***\n");
            Msg_AsFem_Exit();
        }
    }


    string substr=_InputFileName.substr(0,_InputFileName.find_first_of('.'))+"_mesh.vtu";
    _MeshFileName=substr;
    mesh.SetMeshFileName(substr);

    if(!HasType){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type= found in [mesh] block              !!!            ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        type=asfem[gmsh] should be given            !!!            ***\n");   
        Msg_AsFem_Exit();
    }

    
    PetscPrintf(PETSC_COMM_WORLD,"***   start to create mesh ...                                        ***\n");     
    if(!mesh.CreateMesh()){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: create mesh failed !!! please check your input file  !!!   ***\n");
        Msg_AsFem_Exit();
    }
    if(IsSaveMesh){
        mesh.SaveMesh();
        PetscPrintf(PETSC_COMM_WORLD,"***     save mesh to [%45s]! ***\n",substr.c_str());    
    }
    PetscPrintf(PETSC_COMM_WORLD,"***   mesh generation finished !                                      ***\n");     
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