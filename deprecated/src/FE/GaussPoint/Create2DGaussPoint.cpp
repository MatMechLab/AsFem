#include "FE/QPoint.h"

void QPoint::Create2DGaussPoint(MeshType meshtype)
{
    switch (meshtype)
    {
        case MeshType::QUAD4:
        case MeshType::QUAD8:
        case MeshType::QUAD9:
        {
            QPoint qp1d(1,GetQpOrder());
            qp1d.Create1DGaussPoint();
            _nQpPoints=qp1d.GetQpPointsNum()*qp1d.GetQpPointsNum();
            _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);
            int k=0;
            // for this regular type mesh , it is just simple tensor product
            for(int j=0;j<qp1d.GetQpPointsNum();j++)
            {
                for(int i=0;i<qp1d.GetQpPointsNum();i++)
                {
                    k+=1;
                    (*this)(k,1)=qp1d(i+1,1);
                    (*this)(k,2)=qp1d(j+1,1);
                    (*this)(k,0)=qp1d(i+1,0)*qp1d(j+1,0);
                }
            }
            return;
        }
        case MeshType::TRI3:
        case MeshType::TRI6:
        {
            switch (GetQpOrder())
            {
                case 0:
                case 1:
                {
                    _nQpPoints=1;
                    _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);
                    (*this)(1,1)=1.0/3.0;
                    (*this)(1,2)=1.0/3.0;
                    (*this)(1,0)=0.5;
                    return;
                }
                case 2:
                {
                    _nQpPoints=3;
                    _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);
                    (*this)(1,1)=2.0/3.0;
                    (*this)(1,2)=1.0/6.0;
                    (*this)(1,0)=1.0/6.0;

                    (*this)(2,1)=1.0/6.0;
                    (*this)(2,2)=2.0/3.0;
                    (*this)(2,0)=1.0/6.0;
                
                    (*this)(3,1)=1.0/6.0;
                    (*this)(3,2)=1.0/6.0;
                    (*this)(3,0)=1.0/6.0;
                    return;
                }
                case 3:
                {
                    _nQpPoints=4;
                    _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);

                    (*this)(1,1)=double(1.5505102572168219018027159252941e-01L);
                    (*this)(1,2)=double(1.7855872826361642311703513337422e-01L);
                    (*this)(1,0)=double(1.5902069087198858469718450103758e-01L);

                    (*this)(2,1)=double(6.4494897427831780981972840747059e-01L);
                    (*this)(2,2)=double(7.5031110222608118177475598324603e-02L);
                    (*this)(2,0)=double(9.0979309128011415302815498962418e-02L);

                    (*this)(3,1)=double(1.5505102572168219018027159252941e-01L);
                    (*this)(3,2)=double(6.6639024601470138670269327409637e-01L);
                    (*this)(3,0)=double(1.5902069087198858469718450103758e-01L);

                    (*this)(4,1)=double(6.4494897427831780981972840747059e-01L);
                    (*this)(4,2)=double(2.8001991549907407200279599420481e-01L);
                    (*this)(4,0)=double(9.0979309128011415302815498962418e-02L);

                    return;
                }
                default:
                {
                    Msg_GaussPoint_Invalid2DGaussOrder(GetQpOrder());
                    Msg_AsFem_Exit();
                    break;
                }
            }
            return;
        }
        default:
        {
            Msg_GaussPoint_Invalid2DMeshType();
            Msg_AsFem_Exit();
            break;
        }
    }
}

//******************************************************
void QPoint::Create2DGaussLobattoPoint(MeshType meshtype)
{
    switch (meshtype)
    {
        case MeshType::QUAD4:
        case MeshType::QUAD8:
        case MeshType::QUAD9:
        {
            QPoint qp1d(1,GetQpOrder());
            qp1d.SetQPointType("gausslobatto");
            qp1d.Create1DGaussPoint();
            _nQpPoints=qp1d.GetQpPointsNum()*qp1d.GetQpPointsNum();
            _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);
            int k=0;
            // for this regular type mesh , it is just simple tensor product
            for(int i=0;i<qp1d.GetQpPointsNum();i++)
            {
                for(int j=0;j<qp1d.GetQpPointsNum();j++)
                {
                    k+=1;
                    (*this)(k,1)=qp1d(i+1,1);
                    (*this)(k,2)=qp1d(j+1,1);
                    (*this)(k,0)=qp1d(i+1,0)*qp1d(j+1,0);
                }
            }
            return;
        }
        default:
            Msg_GaussPoint_Invalid2DLobattoMeshType();
            Msg_AsFem_Exit();
            break;
    }
}