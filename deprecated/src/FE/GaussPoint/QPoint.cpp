#include "FE/QPoint.h"

QPoint::QPoint()
{
    _qp_coords.clear();
    _nQpPoints=0;
    _nOrder=0;
    _nDim=0;
    _QPointType="gauss";
    _HasSettings=false;
    _HasDim=false;
    _HasOrder=false;
}
//**********************+
QPoint::QPoint(int dim,int order)
{
    _qp_coords.clear();
    _nQpPoints=0;
    _QPointType="gauss";
    _HasSettings=false;
    _HasDim=false;
    _HasOrder=false;

    _nDim=dim;_nOrder=order;
}
