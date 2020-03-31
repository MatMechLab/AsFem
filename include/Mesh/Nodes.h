//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_NODES_H
#define ASFEM_NODES_H

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

class Nodes
{
public:
    Nodes() {_coords.clear();_nNodes=0;}
    Nodes(const int &n) {_coords.resize(4*n,0.0);_nNodes=n;}
    inline void InitNodes(const int &n) {_coords.resize(4*n,0.0);_nNodes=n;}

    inline int GetNodesNum() const {return _nNodes;}
    inline double  operator()(const int &i,const int &j) const {return _coords[(i-1)*4+j];}
    inline double& operator()(const int &i,const int &j) {return _coords[(i-1)*4+j];}

    // operator overload
    inline Nodes operator+(const double &val) const{
        Nodes temp(GetNodesNum());
        for(unsigned int i=0;i<_coords.size();++i) temp._coords[i]=_coords[i]+val;
        return temp;
    }
    inline Nodes& operator+=(const double &val) {
        for(unsigned int i=0;i<_coords.size();++i) _coords[i]=_coords[i]+val;
        return *this;
    }
    inline Nodes operator+(const Nodes &node) const {
        Nodes temp(GetNodesNum());
        for(unsigned int i=0;i<_coords.size();++i) temp._coords[i]=_coords[i]+node._coords[i];
        return temp;
    }
    inline Nodes& operator+=(const Nodes &node) {
        for(unsigned int i=0;i<_coords.size();++i) _coords[i]=_coords[i]+node._coords[i];
        return *this;
    }
    //*********
    inline Nodes operator-(const double &val) const{
        Nodes temp(GetNodesNum());
        for(unsigned int i=0;i<_coords.size();++i) temp._coords[i]=_coords[i]-val;
        return temp;
    }
    inline Nodes& operator-=(const double &val) {
        for(unsigned int i=0;i<_coords.size();++i) _coords[i]=_coords[i]-val;
        return *this;
    }
private:
    int _nNodes;
    vector<double> _coords;// weight+x,y,z
};

#endif // ASFEM_NODES_H