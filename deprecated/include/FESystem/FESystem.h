#ifndef ASFEM_FESYSTEM_H
#define ASFEM_FESYSTEM_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>


#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "BCs/BCSystem.h"

#include "ElmtSystem/ElmtSystem.h"
#include "MaterialSystem/MaterialSystem.h"
#include "SolutionSystem/SolutionSystem.h"

#include "FE/FE.h"
#include "FE/QPoint.h"
#include "FE/QPBlockInfo.h"




#include "Eigen/Eigen"

using namespace std;

class FESystem
{
public:
    FESystem();
    void InitFESystem(Mesh &mesh,
                      DofHandler &dofHandler,
                      //SolutionSystem &solutionSystem,
                      QPBlockInfo &qpBlockInfo);

    FE& GetFEPtr(){return fe;}

    int GetBulkElmtGPointNums() const{return fe._qp_bulk.GetQpPointsNum();}

public:
    void SetHistNumPerGPoint(const int &n){_nHist=n;}
    void SetProjNumPerNode(const int &n) {_nProj=n;}
    void SetKMatrixFactor(const double &factor){_KMatrixFactor=factor;}
    inline double GetKMatrixFactor() const{return _KMatrixFactor;}
    void ResetMaxAMatrixValue(){_MaxKMatrixValue=1.0e3;}
    void SetMaxAMatrixValue(const double &val) {_MaxKMatrixValue=val;}
    inline double GetMaxAMatrixValue()const {return _MaxKMatrixValue;}
    inline double GetVolume() const {return _Volumes;}
    // for FEM simulation related functions
    void FormFE(const int &isw,const double &t,const double &dt,const double (&ctan)[2],
                Mesh &mesh,DofHandler &dofHandler,
                ElmtSystem &elmtSystem,MaterialSystem &mateSystem,
                const Eigen::VectorXd &U,const Eigen::VectorXd &V,
                Eigen::MatrixXd &Hist,const Eigen::MatrixXd &HistOld,Eigen::MatrixXd &Proj,
                Eigen::SparseMatrix<double,Eigen::RowMajor> &AMATRIX,Eigen::VectorXd &RHS);
    
    
private:
    void Projection(const int &nNodes,Eigen::MatrixXd &Proj,const vector<double> &elProj,ShapeFun &shp,const double &DetJac);
    void AssembleLocalHistToGlobal(const int &elmtid,const int &gpInd,const vector<double> &gpHist,Eigen::MatrixXd &Hist);
    
    void AssembleLocalToGlobal(const int &isw,const int &ndofs,const vector<int> &elDofs,
                               const Eigen::MatrixXd &localK,const Eigen::VectorXd &localR,
                               Eigen::SparseMatrix<double,Eigen::RowMajor> &AMATRIX,Eigen::VectorXd &RHS);
    
    void AssembleLocalToGlobal(const int &isw,const int &ndofs,
                               const vector<int> &elDofs,const vector<double> &elDofsActiveFlag,
                               const Eigen::MatrixXd &localK,const Eigen::VectorXd &localR,
                               Eigen::SparseMatrix<double,Eigen::RowMajor> &AMATRIX,Eigen::VectorXd &RHS);

public:
    void PrintFESystemInfo() const;


private:
    FE fe;


private:
    double _Volumes=0.0;
    Eigen::MatrixXd _localK;
    Eigen::VectorXd _localR;
    
    Nodes _elNodes;
    vector<int> _elConn,_elDofs;
    vector<double> _elDofsActiveFlag;
    vector<double> _elU,_elV;
    vector<double> _gpU,_gpV;
    vector<double> _gpHist,_gpHistOld,_gpProj;
    vector<Eigen::Vector3d> _gpGradU,_gpGradV;
    vector<double> _MaterialValues;
    int _nHist,_nProj;
    double _MaxKMatrixValue=1.0e3;
    double _KMatrixFactor=1.0e16;
};

#endif // ASFEM_FESYSTEM_H