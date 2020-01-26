#include "FE/ShapeFun.h"

void ShapeFun::Calc(const double &xi,const Nodes &nodes)
{
    Compute1DLagrangeShapeFun(xi,nodes);
}

//*************************************
void ShapeFun::Calc(const double &xi,const double &eta,const Nodes &nodes)
{
    Compute2DLagrangeShapeFun(xi,eta,nodes);
}
void ShapeFun::Calc(const double &xi,const double &eta,const double &zeta,const Nodes &nodes)
{
    Compute3DLagrangeShapeFun(xi,eta,zeta,nodes);
}