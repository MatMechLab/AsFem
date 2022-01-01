//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.07.08
//+++ Purpose: Implement the reader for [mesh] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/MeshBlockReader.h"


/**
 * Basic usage of [mesh] block:<br>
 * case 1:<br>
 * [mesh]<br>
 *   type=asfem<br>
 *   dim=2 <br>
 *   xmin=0.0 <br>
 *   xmax=1.0 <br>
 *   ymin=0.0 <br>
 *   ymax=1.0 <br>
 *   nx=20 <br>
 *   ny=20 <br>
 *   meshtype=quad4 <br>
 * [end] <br>
 * case 2: <br>
 * [mesh] <br>
 *   type=gmsh <br>
 *   file=mymesh.msh <br>
 *   savemesh=true <br>
 * [end] <br>
 */
void MeshBlockReader::PrintHelper(){
    MessagePrinter::PrintNormalTxt("The complete information for [mesh] block:",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("type=asfem,gmsh2,gmsh4,abaqus,iga",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("dim=1,2,3",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("xmin=0.0",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("xmax=1.0",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("ymin=0.0",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("ymax=0.0",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("zmin=0.0",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("zmax=1.0",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("nx=10",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("ny=10",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("nz=10",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("meshtype=edge2,edge3,quad4,quad8,quad9,hex8,hex20,hex27",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("savemesh=true,false",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("file=meshfile.msh,meshfile.inp",MessageColor::BLUE);

}
//*****************************************************************

/**
 * The the mesh block information as well as the mesh file
 * @param in the file read stream
 * @param str the string which contains '[mesh]' line
 * @param linenum the line number of the '[mesh]' line
 * @param mesh the reference of the mesh class
 */
bool MeshBlockReader::ReadMeshBlock(ifstream &in, string str, int &linenum, Mesh &mesh){
    double xmin=0.0,xmax=1.0;
    double ymin=0.0,ymax=1.0;
    double zmin=0.0,zmax=1.0;
    int dim=2;
    int nx=2,ny=2,nz=2;
    bool IsSuccess=false;
    vector<double> numbers; /**< we use numbers to store the float numbers split from strings */
    bool HasXmin=false,HasXmax=false;
    bool HasYmin=false,HasYmax=false;
    bool HasZmin=false,HasZmax=false;
    bool HasNx=false,HasNy=false,HasNz=false;
    bool IsBuiltIn=false;
    string meshtypename;
    bool HasType=false,HasDim=false,HasMeshType=false;
    bool IsPrint=false,IsDepPrint=false;
    char buff[55]; /**< buff is used to convert chars to strings for the message output*/
    string meshtype;    
    MeshIO meshio; /**< we use the meshio class to read the mesh file and mesh information for our mesh class*/
    bool IsSaveMesh=false;
    // now the str already contains '[mesh]'
    getline(in,str);linenum+=1;
    str=StringUtils::RemoveStrSpace(str);
    meshtypename.clear();

    while(str.find("[end]")==string::npos){
        if(StringUtils::IsCommentLine(str)||str.length()<1){
            // if the current line is a comment line or an empty line, then we skip it 
            getline(in,str);linenum+=1;
            str=StringUtils::RemoveStrSpace(str);
            continue;
        }
        if(StringUtils::StrToLower(str).find("type=helper")!=string::npos){
            PrintHelper();
            return false;
        }
        else if(str.find("type=asfem")!=string::npos){
            HasType=true;
            IsBuiltIn=true;
            getline(in,str);linenum+=1;
            str=StringUtils::StrToLower(str);
            // now we start to read the [mesh] information for the built-in mesh generation
            while(str.find("[end]")==string::npos){
                str=StringUtils::RemoveStrSpace(str);
                str=StringUtils::StrToLower(str);
                if(StringUtils::IsCommentLine(str)||str.length()<1){
                    // if the current line is a comment line or an empty line, then we skip it 
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
                else if(str.find("xmin=")!=string::npos){
                    if(!HasDim){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("xmin= should be given after dim= in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
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
                else if(str.find("xmax=")!=string::npos){
                    if(!HasDim){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("xmax= should be given after dim= in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
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
                else if(str.find("ymin=")!=string::npos){
                    if(!HasDim){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("ymin= should be given after dim= in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
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
                else if(str.find("ymax=")!=string::npos){
                    if(!HasDim){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("ymax= should be given after dim= in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
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
                else if(str.find("zmin=")!=string::npos){
                    if(!HasDim){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("zmin= should be given after dim= in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
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
                else if(str.find("zmax=")!=string::npos){
                    if(!HasDim){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("zmax= should be given after dim= in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
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
                else if(str.find("nx=")!=string::npos){
                    if(!HasDim){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("nx= should be given after dim= in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
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
                else if(str.find("ny=")!=string::npos){
                    if(!HasDim){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("ny= should be given after dim= in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
                    if(dim<2){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("ny= should be given for dim>=2 case in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
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
                            MessagePrinter::PrintErrorTxt("ny= is invalid in [mesh] block, nx= integer value should be given in the [mesh] block");
                            MessagePrinter::AsFem_Exit();
                        }
                        HasNy=true;
                    }
                }
                else if(str.find("nz=")!=string::npos){
                    if(!HasDim){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("nz= should be given after dim= in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
                    if(dim<3){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("nz= should be given for dim=3 case in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
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
                else if(str.find("meshtype=")!=string::npos){
                    if(str.length()<13){
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("invalid meshtyp= information in the [mesh] block");
                        MessagePrinter::AsFem_Exit();
                    }
                    else{
                        meshtypename=str.substr(9,str.length());
                        HasMeshType=true;
                    }
                }
                else if(str.find("savemesh=")!=string::npos){
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
                else if(str.find("printmesh=")!=string::npos){
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
                    MessagePrinter::PrintErrorTxt("if you dont know what to do, please try 'type=helper' ");
                    MessagePrinter::AsFem_Exit();
                }
                getline(in,str);linenum+=1;
            }/** < end of the while loop for 'type=asfem' case */

            // now, since all the information is read for 'type=asfem' case, we try to set up the mesh class
            // and generate the related mesh!
            if(!HasDim){
                MessagePrinter::PrintErrorTxt("dim=1[2,3] should be given for 'type=asfem' case in the [mesh] block");
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
                    MessagePrinter::PrintErrorTxt("for dim=1 case, you only need nx in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(HasYmin||HasYmax||HasZmin||HasZmax){
                    MessagePrinter::PrintErrorTxt("for dim=1 case, you only need xmin,xmax in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasMeshType){
                    MessagePrinter::PrintErrorTxt("meshtype is not assigned, please check your input file");
                    MessagePrinter::AsFem_Exit();
                }

                if(meshtypename!="edge2"&&meshtypename!="edge3"&&meshtypename!="edge4"){
                    MessagePrinter::PrintErrorTxt("unsupported 1d mesh type in the [mesh] block, meshtype=edge2,edge3,edge4 is expected");
                    MessagePrinter::AsFem_Exit();
                }
                IsSuccess=true;
                mesh.SetBulkMeshDim(dim);
                mesh.SetBulkMeshNx(nx);
                mesh.SetBulkMeshXmin(xmin);
                mesh.SetBulkMeshXmax(xmax);
                mesh.SetBulkMeshMeshTypeName(meshtypename);
            }
            else if(dim==2){
                if(!HasNx){
                    MessagePrinter::PrintErrorTxt("nx should be given for dim=2 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasNy){
                    MessagePrinter::PrintErrorTxt("ny should be given for dim=2 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasXmin&&xmin!=0.0){
                    MessagePrinter::PrintErrorTxt("xmin=val should be given for dim=1 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasXmax&&xmax!=1.0){
                    MessagePrinter::PrintErrorTxt("xmax=val should be given for dim=1 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasYmin&&ymin!=0.0) {
                    MessagePrinter::PrintErrorTxt("ymin=val should be given for dim=1 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasYmax&&ymax!=1.0){
                    MessagePrinter::PrintErrorTxt("ymax=val should be given for dim=2 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(HasNz){
                    MessagePrinter::PrintErrorTxt("for dim=2 case, you only need nx and ny in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(HasZmin||HasZmax){
                    MessagePrinter::PrintErrorTxt("for dim=2 case, you only need xmin/xmax and ymin/ymax in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasMeshType){
                    MessagePrinter::PrintErrorTxt("meshtype is not assigned yet, please check your input file");
                    MessagePrinter::AsFem_Exit();
                }
                if(meshtypename!="quad4"&&meshtypename!="quad8"&&meshtypename!="quad9"){
                    MessagePrinter::PrintErrorTxt("unsupported 2d mesh type in the [mesh] block, meshtype=quad4,quad8,quad9 is expected");
                    MessagePrinter::AsFem_Exit();
                }
                IsSuccess=true;
                mesh.SetBulkMeshDim(dim);
                mesh.SetBulkMeshNx(nx);
                mesh.SetBulkMeshNy(ny);
                mesh.SetBulkMeshXmin(xmin);
                mesh.SetBulkMeshXmax(xmax);
                mesh.SetBulkMeshYmin(ymin);
                mesh.SetBulkMeshYmax(ymax);
                mesh.SetBulkMeshMeshTypeName(meshtypename);
            }
            else if(dim==3){
                if(!HasNx){
                    MessagePrinter::PrintErrorTxt("nx should be given for dim=3 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasNy){
                    MessagePrinter::PrintErrorTxt("ny should be given for dim=3 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasNz){
                    MessagePrinter::PrintErrorTxt("nz should be given for dim=3 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasXmin&&xmin!=0.0){
                    MessagePrinter::PrintErrorTxt("xmin=val should be given for dim=3 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasXmax&&xmax!=1.0){
                    MessagePrinter::PrintErrorTxt("xmax=val should be given for dim=3 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasYmin&&ymin!=0.0){
                    MessagePrinter::PrintErrorTxt("ymin=val should be given for dim=3 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasYmax&&ymax!=1.0){
                    MessagePrinter::PrintErrorTxt("ymax=val should be given for dim=3 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasZmin&&zmin!=0.0){
                    MessagePrinter::PrintErrorTxt("zmin=val should be given for dim=3 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }
                if(!HasZmax&&zmax!=1.0){
                    MessagePrinter::PrintErrorTxt("zmax=val should be given for dim=3 case in the [mesh] block");
                    MessagePrinter::AsFem_Exit();
                }

                if(!HasMeshType){
                    MessagePrinter::PrintErrorTxt("meshtype is not assigned yet, please check your input file");
                    MessagePrinter::AsFem_Exit();
                }

                if(meshtypename!="hex8"&&meshtypename!="hex20"&&meshtypename!="hex27"){
                    MessagePrinter::PrintErrorTxt("unsupported 3d mesh type in the [mesh] block, meshtype=hex8,hex20,hex27 is expected");
                    MessagePrinter::AsFem_Exit();
                }
                IsSuccess=true;
                mesh.SetBulkMeshDim(dim);
                mesh.SetBulkMeshNx(nx);
                mesh.SetBulkMeshNy(ny);
                mesh.SetBulkMeshNz(nz);
                mesh.SetBulkMeshXmin(xmin);
                mesh.SetBulkMeshXmax(xmax);
                mesh.SetBulkMeshYmin(ymin);
                mesh.SetBulkMeshYmax(ymax);
                mesh.SetBulkMeshZmin(zmin);
                mesh.SetBulkMeshZmax(zmax);
                mesh.SetBulkMeshMeshTypeName(meshtypename);
            }
            IsBuiltIn=true;
            IsSuccess=false;
            break;
        }
        else if(str.find("type=gmsh")!=string::npos){
            string str0;
            HasType=true;
            IsBuiltIn=false;
            bool HasFileName=false;
            while(str.find("[end]")==string::npos){
                getline(in,str);linenum+=1;
                str=StringUtils::RemoveStrSpace(str);
                str0=str; /** < we always keep str0 as the original string */
                str=StringUtils::StrToLower(str);
                if(StringUtils::IsCommentLine(str)||str.length()<1){
                    // if this is a comment line or empty line, we skip it
                    continue;
                }

                if(str.find("file=")!=string::npos){
                    if(str.compare(str.length()-4,4,".msh")==0){
                        string filename=str0.substr(5,string::npos);
                        meshio.SetMeshFileName(filename);
                        HasFileName=true;
                        MessagePrinter::PrintNormalTxt("Start to import mesh from gmsh ...");
                        IsSuccess=meshio.ReadMeshFromFile(mesh);
                        MessagePrinter::PrintNormalTxt("Import mesh finished !");
                        // mesh.SetMeshMode(false);
                    }
                    else{
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("can not find msh mesh file in your 'file=' option of the [mesh] block, please check your input file");
                        MessagePrinter::AsFem_Exit();
                    }
                }
                else if(str.find("printmesh=")!=string::npos){
                    int i=str0.find_first_of('=');
                    string substr=str0.substr(i+1,str0.length());
                    substr=StringUtils::RemoveStrSpace(substr);
                    if(substr.find("true")!=string::npos){
                        IsPrint=true;IsDepPrint=false;
                    }
                    else if(substr.find("false")!=string::npos){
                        IsPrint=false;IsDepPrint=false;
                    }
                    else if(substr.find("dep")!=string::npos){
                        IsPrint=true;IsDepPrint=true;
                    }
                    else{
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("unsupported option in printmesh= in the [mesh] block, printmesh=true[false,dep] is expected");
                        MessagePrinter::AsFem_Exit();
                    }
                }
                else if(str.find("savemesh=true")!=string::npos){
                    int i=str.find_first_of('=');
                    string substr=str.substr(i+1,str.length());
                    substr=StringUtils::RemoveStrSpace(substr);
                    if(substr.find("true")!=string::npos){
                        IsSaveMesh=true;
                    }
                    else if(substr.find("false")!=string::npos){
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
                else if(str.find("[end]")!=string::npos){
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
                } /**< end of option selection */
                getline(in,str);linenum+=1;
            }/**< end of the while loop */
            if(!HasFileName){
                IsSuccess=false;
                snprintf(buff,55,"line-%d has some errors",linenum);
                MessagePrinter::PrintErrorTxt(string(buff));
                MessagePrinter::PrintErrorTxt("file=correct file name should be given in the [mesh] block");
                MessagePrinter::AsFem_Exit();
            } 
            break;
        }
        else if(str.find("type=abaqus")!=string::npos){
            string str0;
            HasType=true;
            IsBuiltIn=false;
            bool HasFileName=false;
            while(str.find("[end]")==string::npos){
                getline(in,str);linenum+=1;
                str=StringUtils::RemoveStrSpace(str);
                str0=str;
                str=StringUtils::StrToLower(str);
                if(StringUtils::IsCommentLine(str)||str.length()<1){
                    continue;
                }

                if(str.find("file=")!=string::npos){
                    if(str.compare(str.length()-4,4,".inp")==0){
                        string filename=str0.substr(5,string::npos);
                        meshio.SetMeshFileName(filename);
                        MessagePrinter::PrintNormalTxt("Start to import mesh from abaqus ...");
                        IsSuccess=meshio.ReadMeshFromFile(mesh);
                        MessagePrinter::PrintNormalTxt("Import mesh finished !");
                        HasFileName=true;
                    }
                    else{
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("can not find the inp file in the [mesh] block, please check your input file");
                        MessagePrinter::AsFem_Exit();
                    }
                }
                else if(str.find("savemesh=")!=string::npos){
                    int i=str.find_first_of('=');
                    string substr=str.substr(i+1,str.length());
                    substr=StringUtils::RemoveStrSpace(substr);
                    if(substr.find("true")!=string::npos){
                        IsSaveMesh=true;
                    }
                    else if(substr.find("false")!=string::npos){
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
                else if(str.find("printmesh=")!=string::npos){
                    int i=str.find_first_of('=');
                    string substr=str.substr(i+1,str.length());
                    substr=StringUtils::RemoveStrSpace(substr);
                    if(substr.find("true")!=string::npos){
                        IsPrint=true;IsDepPrint=false;
                    }
                    else if(substr.find("false")!=string::npos){
                        IsPrint=false;IsDepPrint=false;
                    }
                    else if(substr.find("dep")!=string::npos){
                        IsPrint=true;IsDepPrint=true;
                    }
                    else{
                        snprintf(buff,55,"line-%d has some errors",linenum);
                        MessagePrinter::PrintErrorTxt(string(buff));
                        MessagePrinter::PrintErrorTxt("unsupported option in printmesh= in the [mesh] block, printmesh=true[false] is expected");
                        MessagePrinter::AsFem_Exit();
                    }
                }
                else if(str.find("[end]")!=string::npos){
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
            break;
        }
        else if(str.find("[end]")!=string::npos){
            break;
        }
        else{
            snprintf(buff,55,"line-%d has some errors",linenum);
            MessagePrinter::PrintErrorTxt(string(buff));
            MessagePrinter::PrintErrorTxt("unknown option in the [mesh] block");
            MessagePrinter::AsFem_Exit();
        }
    }

    if(!HasType){
        MessagePrinter::PrintErrorTxt("no type= found in the [mesh] block, type=asfem[gmsh] should be given");
        MessagePrinter::AsFem_Exit();
    }

    if(IsBuiltIn){
        MessagePrinter::PrintNormalTxt("Start to create mesh ...");
        if(!mesh.CreateMesh()){
            MessagePrinter::PrintErrorTxt("create mesh failed !!! please check your input file");
            MessagePrinter::AsFem_Exit();
        }
        IsSuccess=true;
        MessagePrinter::PrintNormalTxt("Mesh generation finished !");
    }
    if(IsSaveMesh){
        string substr=_InputFileName.substr(0,_InputFileName.find_first_of('.'))+"_mesh.vtu";
        _MeshFileName=substr;
        mesh.SaveMesh(_MeshFileName);
        snprintf(buff,55,"Save mesh to [%39s]",substr.c_str());
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

    MessagePrinter::PrintDashLine();

    return IsSuccess;
}
