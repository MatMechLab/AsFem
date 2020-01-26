#ifndef ASFEM_QPOINT_H
#define ASFEM_QPOINT_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>


//******************************************
#include "MessagePrint/MessagePrint.h"
#include "Mesh/MeshTypeDefine.h"
using namespace std;

class QPoint
{
public:
    QPoint();
    QPoint(int dim,int order);


    void SetDim(int dim) {this->_nDim=dim;}
    void SetQPointOrder(int order) {this->_nOrder=order;}
    void SetQPointType(string type="gauss") {this->_QPointType=type;}
    void CreateQPoints(MeshType meshtype);

    // get function
    inline int GetDim() const {return this->_nDim;}
    inline int GetQpOrder() const {return this->_nOrder;}
    inline int GetQpPointsNum() const {return this->_nQpPoints;}
    
    inline double operator()(int i,int j) const {return this->_qp_coords[(i-1)*(_nDim+1)+j];}
    inline double& operator()(int i,int j) {return this->_qp_coords[(i-1)*(_nDim+1)+j];}

    inline double GetIthQpPointJthCoord(int i,int j) const {return this->_qp_coords[(i-1)*(_nDim+1)+j];}

private:
    void Create1DGaussPoint();
    void Create1DGaussLobattoPoint();

    void Create2DGaussPoint(MeshType meshtype);
    void Create2DGaussLobattoPoint(MeshType meshtype);

    void Create3DGaussPoint(MeshType meshtype);
    void Create3DGaussLobattoPoint(MeshType meshtype);

private:
    vector<double> _qp_coords;
    int _nQpPoints,_nOrder;
    int _nDim;
    string _QPointType;
    bool _HasSettings=false;
    bool _HasDim,_HasOrder;

};

#endif 