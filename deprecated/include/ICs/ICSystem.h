#ifndef ASFEM_ICSYSTEM_H
#define ASFEM_ICSYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "MessagePrint/MessagePrint.h"

#include "ICBlock.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"

#include "Eigen/Eigen"

using namespace std;

class ICSystem
{
public:
    ICSystem();

    void ApplyIC(Mesh &mesh,DofHandler &dofHandler,Eigen::VectorXd &U);

    // setting functions
    bool AddICBlock(const ICBlock &icblock);

    // get information of BCSystem
    inline int GetICBlockNum() const {return _nICBlocks;}
    inline ICBlock GetIthICBlock(int i) {return _ICBlockList[i-1];}


    void PrintICBlockInfo() const;

private:
    void ApplyConstIC(const vector<double> &Params,const int &DofIndex,Mesh &mesh,DofHandler &dofHandler,Eigen::VectorXd &U);
    void ApplyRandIC(const vector<double> &Params,const int &DofIndex,Mesh &mesh,DofHandler &dofHandler,Eigen::VectorXd &U);
    void ApplyCircleIC(const vector<double> &Params,const int &DofIndex,Mesh &mesh,DofHandler &dofHandler,Eigen::VectorXd &U);

private:
    int _nICBlocks;
    vector<ICBlock> _ICBlockList;
};

#endif //