#ifndef ASFEM_ELMTSYSTEM_H
#define ASFEM_ELMTSYSTEM_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "ElmtBlock.h"

#include "FE/ShapeFun.h"

#include "Utils/MathUtils.h"
#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"

#include "ElmtType.h"

using namespace std;

class ElmtSystem
{
public:
    ElmtSystem();

    // setting functions
    bool AddElmtBlock(const ElmtBlock &elmtBlock);

    // get information of BCSystem
    inline int GetElmtBlockNum() const {return _nElmtBlock;}
    inline ElmtBlock GetIthElmtBlock(int i) {return _ElmtBlockList[i-1];}

    void SetIthElmtBlockMateID(const int i,MaterialType mateid){_ElmtBlockList[i-1]._MateTypeID=mateid;}

    void RunElmtLib(const int &isw,const ElmtType &elmttype,
                    const int &nDim,
                    const int &nNodes,
                    const double &JxW,
                    const double &t,const double &dt,const double (&ctan)[2],
                    const Eigen::Vector3d &gpCoord,
                    const vector<double> &gpU,const vector<double> &gpV,
                    const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                    const ShapeFun &shp,
                    const vector<double> &MaterialValues,
                    const vector<RankTwoTensor> &Rank2MaterialValues,
                    const vector<RankFourTensor> &Rank4MaterialValues,
                    vector<double> &Hist,const vector<double> &HistOld,vector<double> &Proj,
                    Eigen::MatrixXd &K,Eigen::VectorXd &rhs);

    
    void PrintElmtBlockInfo() const;

private:
    
    void Poisson(const int &isw,
                const int &nDim,
                const int &nNodes,
                const double &JxW,
                const double &t,const double &dt,const double (&ctan)[2],
                const Eigen::Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &MaterialValues,
                const vector<RankTwoTensor> &Rank2MaterialValues,
                const vector<RankFourTensor> &Rank4MaterialValues,
                vector<double> &Hist,const vector<double> &HistOld,vector<double> &Proj,
                Eigen::MatrixXd &K,Eigen::VectorXd &rhs);

    void Diffusion(const int &isw,
                const int &nDim,
                const int &nNodes,
                const double &JxW,
                const double &t,const double &dt,const double (&ctan)[2],
                const Eigen::Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &MaterialValues,
                const vector<RankTwoTensor> &Rank2MaterialValues,
                const vector<RankFourTensor> &Rank4MaterialValues,
                vector<double> &Hist,const vector<double> &HistOld,vector<double> &Proj,
                Eigen::MatrixXd &K,Eigen::VectorXd &rhs);
    

    void CahnHilliard(const int &isw,
                const int &nDim,
                const int &nNodes,
                const double &JxW,
                const double &t,const double &dt,const double (&ctan)[2],
                const Eigen::Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &MaterialValues,
                const vector<RankTwoTensor> &Rank2MaterialValues,
                const vector<RankFourTensor> &Rank4MaterialValues,
                vector<double> &Hist,const vector<double> &HistOld,vector<double> &Proj,
                Eigen::MatrixXd &K,Eigen::VectorXd &rhs);
        
    void Wave(const int &isw,
                const int &nDim,
                const int &nNodes,
                const double &JxW,
                const double &t,const double &dt,const double (&ctan)[2],
                const Eigen::Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &MaterialValues,
                const vector<RankTwoTensor> &Rank2MaterialValues,
                const vector<RankFourTensor> &Rank4MaterialValues,
                vector<double> &Hist,const vector<double> &HistOld,vector<double> &Proj,
                Eigen::MatrixXd &K,Eigen::VectorXd &rhs);

    void Mechanics(const int &isw,
                const int &nDim,
                const int &nNodes,
                const double &JxW,
                const double &t,const double &dt,const double (&ctan)[2],
                const Eigen::Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &MaterialValues,
                const vector<RankTwoTensor> &Rank2MaterialValues,
                const vector<RankFourTensor> &Rank4MaterialValues,
                vector<double> &Hist,const vector<double> &HistOld,vector<double> &Proj,
                Eigen::MatrixXd &K,Eigen::VectorXd &rhs);
    
    void ThermalMechanics(const int &isw,
                const int &nDim,
                const int &nNodes,
                const double &JxW,
                const double &t,const double &dt,const double (&ctan)[2],
                const Eigen::Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &MaterialValues,
                const vector<RankTwoTensor> &Rank2MaterialValues,
                const vector<RankFourTensor> &Rank4MaterialValues,
                vector<double> &Hist,const vector<double> &HistOld,vector<double> &Proj,
                Eigen::MatrixXd &K,Eigen::VectorXd &rhs);

    //*********************************************************
    //*** For Miehe's phase field fracture model
    void LinearElasticPhaseFieldFracture(const int &isw,
                const int &nDim,
                const int &nNodes,
                const double &JxW,
                const double &t,const double &dt,const double (&ctan)[2],
                const Eigen::Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &MaterialValues,
                const vector<RankTwoTensor> &Rank2MaterialValues,
                const vector<RankFourTensor> &Rank4MaterialValues,
                vector<double> &Hist,const vector<double> &HistOld,vector<double> &Proj,
                Eigen::MatrixXd &K,Eigen::VectorXd &rhs);

    //*********************************************************
    //*** For user defined elements (UEL)
    //*********************************************************
    void UserElmt1(const int &isw,
                const int &nDim,
                const int &nNodes,
                const double &JxW,
                const double &t,const double &dt,const double (&ctan)[2],
                const Eigen::Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &MaterialValues,
                const vector<RankTwoTensor> &Rank2MaterialValues,
                const vector<RankFourTensor> &Rank4MaterialValues,
                vector<double> &Hist,const vector<double> &HistOld,vector<double> &Proj,
                Eigen::MatrixXd &K,Eigen::VectorXd &rhs);


private:
    int _nElmtBlock;
    vector<ElmtBlock> _ElmtBlockList;

    Eigen::Vector3d _grad_test,_grad_phi;
};

#endif 