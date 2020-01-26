#include "FE/QPoint.h"


void QPoint::CreateQPoints(MeshType meshtype)
{
    if(GetDim()==1)
    {
        if(_QPointType=="gauss")
        {
            Create1DGaussPoint();
        }
        else if(_QPointType=="gausslobatto")
        {
            Create1DGaussLobattoPoint();
        }
    }
    else if(GetDim()==2)
    {
        if(_QPointType=="gauss")
        {
            Create2DGaussPoint(meshtype);
        }
        else if(_QPointType=="gausslobatto")
        {
            Create2DGaussLobattoPoint(meshtype);
        }
    }
    else if(GetDim()==3)
    {
        if(_QPointType=="gauss")
        {
            Create3DGaussPoint(meshtype);
        }
        else if(_QPointType=="gausslobatto")
        {
            Create3DGaussLobattoPoint(meshtype);
        }
    }
}