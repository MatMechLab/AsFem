//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FE/QPoint.h"


void QPoint::Create2DGaussPoint(MeshType meshtype){
    switch (meshtype){
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
            for(int j=0;j<qp1d.GetQpPointsNum();j++){
                for(int i=0;i<qp1d.GetQpPointsNum();i++){
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
            // for the details, one is referred to:
            // http://www.ce.memphis.edu/7111/notes/class_notes/chapter_03d_slides.pdf
            switch (GetQpOrder())
            {
                case 0:
                case 1:
                {
                    _nQpPoints=1;
                    _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);
                    (*this)(1,1)=1.0/3.0;
                    (*this)(1,2)=1.0/3.0;
                    (*this)(1,0)=1.0;
                    return;
                }
                case 2:
                {
                    _nQpPoints=3;
                    _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);
                    (*this)(1,1)=2.0/3.0;
                    (*this)(1,2)=1.0/6.0;
                    (*this)(1,0)=1.0/3.0;

                    (*this)(2,1)=1.0/6.0;
                    (*this)(2,2)=2.0/3.0;
                    (*this)(2,0)=1.0/3.0;
                
                    (*this)(3,1)=1.0/6.0;
                    (*this)(3,2)=1.0/6.0;
                    (*this)(3,0)=1.0/3.0;
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
                case 4:
                {
                    _nQpPoints=6;
                    _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);

                    (*this)(1,1)= 0.0915762135;
                    (*this)(1,2)= 0.8168475730;
                    (*this)(1,0)= 0.1099517437;

                    (*this)(2,1)= 0.0915762135;
                    (*this)(2,2)= 0.0915762135;
                    (*this)(2,0)= 0.1099517437;

                    (*this)(3,1)= 0.8168475730;
                    (*this)(3,2)= 0.0915762135;
                    (*this)(3,0)= 0.1099517437;

                    (*this)(4,1)= 0.4459484909;
                    (*this)(4,2)= 0.1081030182;
                    (*this)(4,0)= 0.2233815897;

                    (*this)(5,1)= 0.4459484909;
                    (*this)(5,2)= 0.4459484909;
                    (*this)(5,0)= 0.2233815897;

                    (*this)(6,1)= 0.1081030182;
                    (*this)(6,2)= 0.4459484909;
                    (*this)(6,0)= 0.2233815897;

                    return;
                }
                case 5:
                {
                    _nQpPoints=7;
                    _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);

                    (*this)(1,1)= 1.0/3.0;
                    (*this)(1,2)= 1.0/3.0;
                    (*this)(1,0)= 0.2250000000;

                    (*this)(2,1)= 0.1012865073;
                    (*this)(2,2)= 0.7974269854;
                    (*this)(2,0)= 0.1259391805;

                    (*this)(3,1)= 0.1012865073;
                    (*this)(3,2)= 0.1012865073;
                    (*this)(3,0)= 0.1259391805;

                    (*this)(4,1)= 0.7974269854;
                    (*this)(4,2)= 0.1012865073;
                    (*this)(4,0)= 0.1259391805;

                    (*this)(5,1)= 0.0597158718;
                    (*this)(5,2)= 0.4701420641;
                    (*this)(5,0)= 0.1323941528;

                    (*this)(6,1)= 0.4701420641;
                    (*this)(6,2)= 0.4701420641;
                    (*this)(6,0)= 0.1323941528;

                    (*this)(7,1)= 0.4701420641;
                    (*this)(7,2)= 0.0597158718;
                    (*this)(7,0)= 0.1323941528;

                    return;
                }
                default:
                {
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid gauss integration order(=%3d) for 2D case    !!!   ***\n",GetQpOrder());
                    Msg_AsFem_Exit();
                    break;
                }
            }
            return;
        }
        default:
        {
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type for 2D gauss integration       !!!   ***\n");
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
            for(int i=0;i<qp1d.GetQpPointsNum();i++){
                for(int j=0;j<qp1d.GetQpPointsNum();j++){
                    k+=1;
                    (*this)(k,1)=qp1d(i+1,1);
                    (*this)(k,2)=qp1d(j+1,1);
                    (*this)(k,0)=qp1d(i+1,0)*qp1d(j+1,0);
                }
            }
            return;
        }
        default:
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type for 2D gauss lobatto integration !!! ***\n");
            Msg_AsFem_Exit();
            break;
    }
}