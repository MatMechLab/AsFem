#include "FE/ShapeFun.h"

ShapeFun::ShapeFun()
{
    _nOrder=0;_nFuns=0;
    _nDim=0;
    _MeshType=MeshType::NULLTYPE;
    _IsCartesianDeriv=true;// d/dX,d/dY; false-->d/dxi,d/deta...
    _HasDim=false;_HasOrder=false;_HasMeshType=false;
    _values.clear();
    _nValues=0;

    _shape_value.clear();
    _shape_grad.clear();
}

ShapeFun::ShapeFun(int dim,MeshType meshtype)
{
    _nFuns=0;
    _nDim=dim;
    _MeshType=meshtype;
    _IsCartesianDeriv=true;// d/dX,d/dY; false-->d/dxi,d/deta...
    _values.clear();
    _HasDim=true;_HasOrder=true;_HasMeshType=true;
    _nValues=0;

    _shape_value.clear();
    _shape_grad.clear();
}


void ShapeFun::PreCalc()
{
    if(_HasDim&&_HasMeshType)
    {
        switch (GetMeshType())
        {
            case MeshType::EDGE2:
                _nFuns=2;
                _nOrder=1;
                break;
            case MeshType::EDGE3:
                _nFuns=3;
                _nOrder=2;
                break;
            case MeshType::EDGE4:
                _nFuns=4;
                _nOrder=3;
                break;
            case MeshType::QUAD4:
                _nFuns=4;
                _nOrder=1;
                break;
            case MeshType::QUAD8:
                _nFuns=8;
                _nOrder=2;
                break;
            case MeshType::QUAD9:
                _nFuns=9;
                _nOrder=2;
                break;
            case MeshType::TRI3:
                _nFuns=3;
                _nOrder=1;
                break;
            case MeshType::TRI6:
                _nFuns=6;
                _nOrder=2;
                break;
            case MeshType::HEX8:
                _nFuns=8;
                _nOrder=1;
                break;
            case MeshType::HEX20:
                _nFuns=20;
                _nOrder=2;
                break;
            case MeshType::HEX27:
                _nFuns=27;
                _nOrder=2;
                break;
            case MeshType::TET4:
                _nFuns=4;
                _nOrder=1;
                break;
            case MeshType::TET10:
                _nFuns=10;
                _nOrder=2;
                break;
            default:
                Msg_ShapeFun_UnsupportedMeshType();
                Msg_AsFem_Exit();
                break;
        }
        _nValues=(1+GetDim())*_nFuns;
        _values.reserve(_nValues);
        for(int i=0;i<_nValues;++i) _values.push_back(0.0);


        _shape_value.reserve(_nFuns);
        _shape_grad.reserve(_nFuns);
        for(int i=0;i<_nFuns;++i){
            _shape_value.push_back(0.0);
            _shape_grad.push_back(Eigen::Vector3d::Zero());
        }
        for(int i=0;i<_nFuns;++i){
            _shape_grad[i].setZero();
        }
    }
    else{
        Msg_ShapeFun_InformationNotComlete();
    }
}