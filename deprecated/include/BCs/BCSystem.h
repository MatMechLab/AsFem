#ifndef ASFEM_BCSYSTEM_H
#define ASFEM_BCSYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "BCBlock.h"

#include "Mesh/Mesh.h"
#include "Mesh/Nodes.h"
#include "DofHandler/DofHandler.h"
#include "FE/FE.h"



using namespace std;


class Mesh;
class DofHandler;
class FE;


class BCSystem
{
public:
    BCSystem();
    void Init(Mesh &mesh);
    void SetMaxKMatrixValue(const double &val){_MaxKMatrixValue=val;}

    // setting functions
    bool AddBCBlock(const BCBlock &bcblock);

    // get information of BCSystem
    inline int GetBCBlockNum() const {return _nBCBlocks;}
    inline BCBlock GetIthBCBlock(int i) {return _BCBlockList[i-1];}


    void PrintBCBlockInfo() const;


public:
    void ApplyBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,const double &t,const double (&ctan)[2],
                 Eigen::SparseMatrix<double,Eigen::RowMajor> &K,
                 Eigen::VectorXd &RHS,
                 Eigen::VectorXd &U);
                
    void ApplyPreSetDirichletBC(const double &t,Mesh &mesh,DofHandler &dofHander,Eigen::VectorXd &U);
    // void PresetDirichletBC(Mesh &mesh,DofHandler &dofHander,Eigen::VectorXd &U);
    // void ApplyBC(Mesh &mesh,DofHandler &dofHander,
    //              ShapeFun &shp,);

private:
    void ApplyNeumannBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,
                        const int &DofIndex,const double &bcvalue,string bcname,
                        Eigen::VectorXd &RHS);

    void ApplyDirichletBC(Mesh &mesh,DofHandler &dofHandler,
                        const int &DofIndex,const double &bcvalue,string bcname,
                        Eigen::SparseMatrix<double,Eigen::RowMajor> &K,
                        Eigen::VectorXd &RHS,
                        Eigen::VectorXd &U);
    // void ApplyDynamicDirichletBC()
    void ApplyPresetBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,const double (&ctan)[2],
                        const int &DofIndex,const double &bcvalue,string bcname,
                        Eigen::SparseMatrix<double,Eigen::RowMajor> &K,
                        Eigen::VectorXd &RHS,
                        Eigen::VectorXd &U);

private:
    int _nBCBlocks;
    vector<BCBlock> _BCBlockList;

    int _nNodesPerBCElmt=1,_nDim=1;
    double _xi,_eta,_JxW;
    double _Area,_ElArea;

    int _DofIndex;
    double _bcvalue;
    string _bcname;
    Nodes _elNodes;

    double _MaxKMatrixValue=1.0;

};

#endif // ASFEM_BCSYSTEM_H