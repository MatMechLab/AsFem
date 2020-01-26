#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::RunElmtLib(const int &isw,const ElmtType &elmttype,
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
                    Eigen::MatrixXd &K,Eigen::VectorXd &rhs){
    switch (elmttype)
    {
        case ElmtType::POISSON:
            Poisson(isw,nDim,nNodes,JxW,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                    shp,
                    MaterialValues,Rank2MaterialValues,Rank4MaterialValues,
                    Hist,HistOld,Proj,K,rhs);
            break;
        case ElmtType::DIFFUSION:
            Diffusion(isw,nDim,nNodes,JxW,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                      shp,
                      MaterialValues,Rank2MaterialValues,Rank4MaterialValues,
                      Hist,HistOld,Proj,K,rhs);
            break;
        case ElmtType::WAVE:
            Wave(isw,nDim,nNodes,JxW,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                      shp,
                      MaterialValues,Rank2MaterialValues,Rank4MaterialValues,
                      Hist,HistOld,Proj,K,rhs);
            break;
        case ElmtType::CAHNHILLIARD:
            CahnHilliard(isw,nDim,nNodes,JxW,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                         shp,
                         MaterialValues,Rank2MaterialValues,Rank4MaterialValues,
                         Hist,HistOld,Proj,K,rhs);
            break;
        case ElmtType::MECHANICS:
            Mechanics(isw,nDim,nNodes,JxW,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                      shp,
                      MaterialValues,Rank2MaterialValues,Rank4MaterialValues,
                      Hist,HistOld,Proj,K,rhs);
            break;
        case ElmtType::THERMALMECHANICS:
            ThermalMechanics(isw,nDim,nNodes,JxW,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                    shp,
                    MaterialValues,Rank2MaterialValues,Rank4MaterialValues,
                    Hist,HistOld,Proj,K,rhs);
            break;
        case ElmtType::LINEARELASTICPHASEFIELDFRACTURE:
            LinearElasticPhaseFieldFracture(isw,nDim,nNodes,JxW,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                    shp,
                    MaterialValues,Rank2MaterialValues,Rank4MaterialValues,
                    Hist,HistOld,Proj,K,rhs);
            break;
        case ElmtType::USER1:
            UserElmt1(isw,nDim,nNodes,JxW,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                    shp,
                    MaterialValues,Rank2MaterialValues,Rank4MaterialValues,
                    Hist,HistOld,Proj,K,rhs);
            break;
        default:
            break;
    }
}