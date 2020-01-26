#ifndef ASFEM_FE_H
#define ASFEM_FE_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>


//**********************************
//*** AsFem's own header file    ***
//**********************************
#include "QPoint.h"
#include "ShapeFun.h"
#include "Mesh/Nodes.h"
#include "Mesh/Mesh.h"

#include "QPBlockInfo.h"

class FE
{
public:
    FE();
    void SetDim(int dim) {this->_nDim=dim;}
    void SetDimMin(int dim) {this->_nDimMin=dim;}
    void SetOrder(int order) {this->_nOrder=order;}
    void InitFE(Mesh &mesh,string qpointtype);// init the gauss point system and shape fun system
                            // before real calculcation
                            // 

    inline int GetDim() const {return this->_nDim;}
    inline int GetDimMin() const {return this->_nDimMin;}
    inline int GetOrder() const {return this->_nOrder;}



public:
    bool _IsInit=false;
    int _nDim,_nDimMin,_nOrder;
    QPoint _qp_bulk,_qp_surface,_qp_line;
    ShapeFun _shp_bulk,_shp_surface,_shp_line;
    Nodes _nodes,_surface_nodes,_line_nodes;
    
};

#endif // ASFEM_FE_H