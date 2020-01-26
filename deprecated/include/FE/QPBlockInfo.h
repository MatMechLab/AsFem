#ifndef ASFEM_QPBLOCKINFO_H
#define ASFEM_QPBLOCKINFO_H

#include <iostream>
#include <string>
#include <iomanip>

using namespace std;

class QPBlockInfo{
public:
    int _nQpOrder=1;
    bool _SetFromInput=false;// if not qpblock is taken from input
                           // the default order=mesh's order+1
    string _QpType="gauss";// gauss or gausslobatto

    void Reset(){
        _nQpOrder=1;
        _SetFromInput=false;
        _QpType="gauss";
    }
};


#endif