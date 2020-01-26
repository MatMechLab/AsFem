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

    // now str should contain "[mesh]"
    getline(in,str);linenum+=1;
    str=RemoveStrSpace(str);
    if(str.find("type=asfem")!=string::npos||
       str.find("type=AsFem")!=string::npos||
       str.find("type=ASFEM")!=string::npos||
       str.find("type=Asfem")!=string::npos){
        HasType=true;
        IsBuiltIn=true;
        getline(in,str);linenum+=1;
        // read the format for built-in
        while(str.find("[end]")==string::npos&&
              str.find("[END]")==string::npos){
            str=RemoveStrSpace(str);
            if(IsCommentLine(str) || str.length()<1) {
                getline(in,str);linenum+=1;
                continue;
            }
            if(str.find("dim=")!=string::npos||
               str.find("DIM=")!=string::npos||
               str.find("Dim=")!=string::npos){
                numbers=SplitStrNum(str);
                if(numbers.size()<1){
                    Msg_Input_LineError(linenum);
                    cout<<"*** dim=1 [2,3] should be given !!!                        ***"<<endl;
                    Msg_AsFem_Exit();
                }
                else{
                    dim=int(numbers[0]);
                    if(dim<1||dim>3){
                        Msg_Input_LineError(linenum);
                        cout<<"*** dim=1 [2,3] should be given !!!                        ***"<<endl;
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
                    cout<<"*** xmin = value should be given !!!                       ***"<<endl;
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
                    cout<<"*** xmax = value should be given !!!                       ***"<<endl;
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
                    cout<<"*** ymin = value should be given !!!                       ***"<<endl;
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
                    cout<<"*** ymax = value should be given !!!                       ***"<<endl;
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
                    cout<<"*** zmin = value should be given !!!                       ***"<<endl;
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
                    cout<<"*** zmax = value should be given !!!                       ***"<<endl;
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
                    cout<<"*** nx = integer should be given !!!                       ***"<<endl;
                    Msg_AsFem_Exit();
                }
                else{
                    nx=int(numbers[0]);
                    if(nx<1){
                        Msg_Input_LineError(linenum);
                        cout<<"*** nx = integer should be given !!!                       ***"<<endl;
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
                    cout<<"*** ny = integer should be given !!!                       ***"<<endl;
                    Msg_AsFem_Exit();
                }
                else{
                    ny=int(numbers[0]);
                    if(ny<1){
                        cout<<"*** ny = integer should be given !!!                       ***"<<endl;
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
                    cout<<"*** nz = integer should be given !!!                       ***"<<endl;
                    Msg_AsFem_Exit();
                }
                else{
                    nz=int(numbers[0]);
                    if(nz<1)
                    {
                        cout<<"*** nz = integer should be given !!!                       ***"<<endl;
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
                    cout<<"*** invalid meshtype = information !!!                     ***"<<endl;
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
    
                    cout<<"*** Error: unsupported option in savemesh= in [mesh] !!!   ***"<<endl;
                    Msg_AsFem_Exit();
                }
            }
            else if(str.find("[end]")!=string::npos||str.find("[END]")!=string::npos)
            {
                break;
            }
            else{
                Msg_Input_LineError(linenum);
                cout<<"*** Error: unknown option in [mesh] block     !!!          ***"<<endl;
                Msg_AsFem_Exit();
            }
            getline(in,str);linenum+=1;
        }

        //*****************************************************
        if(!HasDim)
        {
            cout<<"*** Error: dim=1[2,3] should be given in mesh block!!      ***"<<endl;
            Msg_AsFem_Exit();
        }

        if(dim==1)
        {
            if(!HasNx)
            {
                cout<<"*** Error: nx should be given for dim=1 case !!!           ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasXmin&&xmin!=0.0)
            {
                cout<<"*** Error: xmin=val should be given for dim=1 case !!!     ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasXmax&&xmax!=1.0)
            {
                cout<<"*** Error: xmax=val should be given for dim=1 case !!!     ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(HasNy||HasNz)
            {
                cout<<"*** Error: for dim=1 case, you only need nx !!!            ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(HasYmin||HasYmax||HasZmin||HasZmax)
            {
                cout<<"*** Error: for dim=1 case, you only need xmin,xmax!!!      ***"<<endl;
                Msg_AsFem_Exit();
            }

            if(meshtype!="edge2"&&meshtype!="edge3"&&meshtype!="edge4")
            {
                Msg_Input_Invalid1DMeshType();
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
        else if(dim==2)
        {
            if(!HasNx)
            {
                cout<<"*** Error: nx should be given for dim=2 case !!!           ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasNy)
            {
                cout<<"*** Error: ny should be given for dim=2 case !!!           ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasXmin&&xmin!=0.0)
            {
                cout<<"*** Error: xmin=val should be given for dim=1 case !!!     ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasXmax&&xmax!=1.0)
            {
                cout<<"*** Error: xmax=val should be given for dim=1 case !!!     ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasYmin&&ymin!=0.0)
            {
                cout<<"*** Error: ymin=val should be given for dim=2 case !!!     ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasYmax&&ymax!=1.0)
            {
                cout<<"*** Error: ymax=val should be given for dim=2 case !!!     ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(HasNz)
            {
                cout<<"*** Error: for dim=2 case, you only need nx,ny !!!         ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(HasZmin||HasZmax)
            {
                cout<<"*** Error: for dim=2 case, you only need x,ymin/max!!!     ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(meshtype!="quad4"&&meshtype!="quad8"&&meshtype!="quad9")
            {
                Msg_Input_Invalid2DMeshType();
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
        else if(dim==3)
        {
            if(!HasNx)
            {
                cout<<"*** Error: nx should be given for dim=3 case !!!           ***"<<endl;
                 Msg_AsFem_Exit();
            }
            if(!HasNy)
            {
                cout<<"*** Error: ny should be given for dim=3 case !!!           ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasNy)
            {
                cout<<"*** Error: nz should be given for dim=3 case !!!           ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasXmin&&xmin!=0.0)
            {
                cout<<"*** Error: xmin=val should be given for dim=3 case !!!     ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasXmax&&xmax!=1.0)
            {
                cout<<"*** Error: xmax=val should be given for dim=3 case !!!     ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasYmin&&ymin!=0.0)
            {
                cout<<"*** Error: ymin=val should be given for dim=3 case !!!     ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasYmax&&ymax!=1.0)
            {
                cout<<"*** Error: ymax=val should be given for dim=3 case !!!     ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasZmin&&zmin!=0.0)
            {
                cout<<"*** Error: zmin=val should be given for dim=3 case !!!     ***"<<endl;
                Msg_AsFem_Exit();
            }
            if(!HasZmax&&zmax!=1.0)
            {
                cout<<"*** Error: zmax=val should be given for dim=3 case !!!     ***"<<endl;
                Msg_AsFem_Exit();
            }

            if(meshtype!="hex8"&&meshtype!="hex20"&&meshtype!="hex27")
            {
                Msg_Input_Invalid3DMeshType();
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
            str.find("type=GMSH")!=string::npos)
    {
        HasType=true;
        IsBuiltIn=false;
        getline(in,str);linenum+=1;
        bool HasFileName=false;
        while(str.find("[end]")==string::npos&&
              str.find("[END]")==string::npos)
        {
            str=RemoveStrSpace(str);
            if(IsCommentLine(str)||str.length()<1)
            {
                getline(in,str);linenum+=1;
                continue;
            }
            if(str.find("file=")!=string::npos||
               str.find("FILE=")!=string::npos)
            {
                if(str.compare(str.length()-4,4,".msh")==0||
                   str.compare(str.length()-4,4,".Msh")==0||
                   str.compare(str.length()-4,4,".MSH")==0)
                {
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
    
                    cout<<"*** Error: unsupported option in savemesh= in [mesh] !!!   ***"<<endl;
                    Msg_AsFem_Exit();
                }
            }
            else if(str.find("[end]")!=string::npos||
                    str.find("[END]")!=string::npos)
            {
                break;
            }
            else if(str.find("[]")!=string::npos)
            {
                Msg_Input_LineError(linenum);
                Msg_Input_BlockBracketNotComplete();
                Msg_AsFem_Exit();
            }
            else{
                Msg_Input_LineError(linenum);
                cout<<"*** Error: unknown option in [mesh] block       !!!  ***"<<endl;
                Msg_AsFem_Exit();
            }
            getline(in,str);linenum+=1;
        }
        if(!HasFileName)
        {
            IsSuccess=false;
            cout<<"*** Error: file=correct file name should be given !!!      ***"<<endl;
            Msg_AsFem_Exit();
        }
    }


    string substr=_InputFileName.substr(0,_InputFileName.find_first_of('.'))+"_mesh.vtu";
    mesh.SetMeshFileName(substr);

    if(!HasType){
        cout<<"*** Error: no type= found in [mesh] block !!!              ***"<<endl;
        cout<<"***        type=asfem[gmsh] should be given !!!            ***"<<endl;
        Msg_AsFem_Exit();
    }

    
    printf("***   start to crate mesh ...                              ***\n");
    mesh.CreateMesh();
    if(IsSaveMesh){
        mesh.SaveMesh();
        printf("***     save mesh to %34s    ***\n",substr.c_str());
    }
    printf("***   mesh generation finished !                           ***\n");

    return IsSuccess;

}