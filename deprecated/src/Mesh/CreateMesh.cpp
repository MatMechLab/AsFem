#include "Mesh/Mesh.h"

void Mesh::CreateMesh()
{
    if(_IsBuiltInMesh){
        if(GetDim()==1){
            Create1DMesh();
        }
        else if(GetDim()==2){
            Create2DMesh();
        }
        else if(GetDim()==3){
            Create3DMesh();
        }
    }
    else{
        if(_UseGmsh){
            ReadMeshFromGmsh();
        }
    }
}