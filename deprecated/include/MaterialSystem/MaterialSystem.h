#ifndef ASFEM_MATERIALSYSTEM_H
#define ASFEM_MATERIALSYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
// #include <map>


#include "MessagePrint/MessagePrint.h"
#include "MaterialBlock.h"
#include "MaterialType.h"

#include "Eigen/Eigen"

#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"

using namespace std;


class MaterialSystem
{
public:
    MaterialSystem();
    void ZeroMateValues(){
        fill(_MateValues.begin(),_MateValues.end(),0.0);
        for(auto &it:_Rank2TensorMateValues) it.SetToZeros();
        for(auto &it:_Rank4TensorMateValues) it.SetToZeros();
    }
    void InitMateValues(){
        _MateValues.resize(_nMateValues,0.0);
        _Rank2TensorMateValues.resize(_nRank2MateValues,RankTwoTensor(0.0));
        _Rank4TensorMateValues.resize(_nRank4MateValues,RankFourTensor(0.0));
    }



    // setting functions
    bool AddMaterialBlock(const MaterialBlock &mateblock);
    void SetMateValNums(const int &nmatevals) {_nMateValues=nmatevals;}
    void SetRank2MateValNums(const int &nmatevals) {_nRank2MateValues=nmatevals;}
    void SetRank4MateValNums(const int &nmatevals) {_nRank4MateValues=nmatevals;}
    // get information of BCSystem
    inline int GetMateValNums() const {return _nMateValues;}
    inline int GetRank2MateValNums() const {return _nRank2MateValues;}
    inline int GetRank4MateValNums() const {return _nRank4MateValues;}
    inline int GetMateBlockNum() const {return _nMaterialBlocks;}
    inline MaterialBlock GetIthMateBlock(const int &i) {return _MaterialBlockList[i-1];}
    inline MaterialType GetIthMateBlockMateTypeID(const int &i)const{
        return _MaterialBlockList[i-1]._MateTypeID;
    }
    void PrintMateBlockInfo() const;


public:
    vector<double> _MateValues;
    vector<RankTwoTensor> _Rank2TensorMateValues;
    vector<RankFourTensor> _Rank4TensorMateValues;

private:
    int _nMaterialBlocks;
    vector<MaterialBlock> _MaterialBlockList;
    int _nMateValues=50;
    int _nRank2MateValues=6;
    int _nRank4MateValues=3;


public:
// vector<int> _elConn,_elDofs;
    // vector<double> _elU,_elV;
    // vector<double> _gpU,_gpV;
    // vector<double> _elHist,_elHistOld,_elProj;
    // vector<Eigen::Vector3d> _gpGradU,_gpGradV;
    void RunMateLib(const MaterialType &matetype,
                    const int &iblock,
                    const int &nDim,
                    const double &t,const double &dt,
                    const Eigen::Vector3d &gpCoord,
                    const vector<double> &gpU,const vector<double> &gpV,
                    const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                    vector<double> &gpHist,const vector<double> &gpHistOld);
private:
    void ConstPoissonMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    void NonlinearPoissonMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    void ConstDiffusionMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    void NonlinearDiffusionMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    //*** for general wave equation's material
    void GenralWaveMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);

    //*** for the free energies and its derivatives material in CH equation
    void CahnHilliardMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);

    //*** for solid mechanics material (stress,strain and constitutive law)
    void LinearElasticMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    // for Miehe's model in linear elastic case
    void LinearElasticPhaseFieldFracMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    void ModifyElasticFracMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);

    //*** Finite strain material
    void FiniteStrainMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    //*** Saint-Venant material
    void SaintVenantMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);
    //*** compressible Neo-Hookean material
    void NeoHookeanMate(const int &nDim,const double &t,const double &dt,
                          const vector<double> InputParams,
                          const Eigen::Vector3d &gpCoord,
                          const vector<double> &gpU,const vector<double> &gpV,
                          const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                          vector<double> &gpHist,const vector<double> &gpHistOld);

    //*** for thermal-mechanics material (stress,strain and constitutive law)
    void ThermalElasticMate(const int &nDim,const double &t,const double &dt,
                        const vector<double> InputParams,
                        const Eigen::Vector3d &gpCoord,
                        const vector<double> &gpU,const vector<double> &gpV,
                        const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                        vector<double> &gpHist,const vector<double> &gpHistOld);

    //*****************************************
    //*** User defined materials(umat)      ***
    //*****************************************
    void UserMate1(const int &nDim,const double &t,const double &dt,
                   const vector<double> InputParams,
                   const Eigen::Vector3d &gpCoord,
                   const vector<double> &gpU,const vector<double> &gpV,
                   const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                   vector<double> &gpHist,const vector<double> &gpHistOld);


                   
};


#endif // ASFEM_MATERIALSYSTEM_H