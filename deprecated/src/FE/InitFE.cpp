#include "FE/FE.h"

void FE::InitFE(Mesh &mesh,string qpointtype)
{
    _shp_bulk=ShapeFun(mesh.GetDim(),mesh.GetBulkElmtType());
    _shp_bulk.PreCalc();

    _qp_bulk=QPoint(mesh.GetDim(),GetOrder());
    _qp_bulk.SetQPointType(qpointtype);
    _qp_bulk.CreateQPoints(mesh.GetBulkElmtType());
    

    _nodes=Nodes(mesh.GetNodesNumPerBulkElmt());
    if(mesh.GetDim()==3){
        _shp_surface=ShapeFun(2,mesh.GetSurfaceElmtType());
        _shp_surface.PreCalc();

        _shp_line=ShapeFun(1,mesh.GetLineElmtType());
        _shp_line.PreCalc();

        _qp_surface=QPoint(2,GetOrder());
        _qp_surface.SetQPointType(qpointtype);
        _qp_surface.CreateQPoints(mesh.GetSurfaceElmtType());

        _qp_line=QPoint(1,GetOrder());
        _qp_line.SetQPointType(qpointtype);
        _qp_line.CreateQPoints(mesh.GetLineElmtType());

        _surface_nodes=Nodes(mesh.GetNodesNumPerSurfaceElmt());
        _line_nodes=Nodes(mesh.GetNodesNumPerLineElmt());
    }
    else if(mesh.GetDim()==2)
    {
        _shp_line=ShapeFun(1,mesh.GetLineElmtType());
        _shp_line.PreCalc();

        _qp_line=QPoint(1,GetOrder());
        _qp_line.SetQPointType(qpointtype);
        _qp_line.CreateQPoints(mesh.GetLineElmtType());

        _line_nodes=Nodes(mesh.GetNodesNumPerLineElmt());
    }
}