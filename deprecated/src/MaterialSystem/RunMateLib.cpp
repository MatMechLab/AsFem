#include "MaterialSystem/MaterialSystem.h"

void MaterialSystem::RunMateLib(const MaterialType &matetype,
                    const int &iblock,
                    const int &nDim,
                    const double &t,const double &dt,
                    const Eigen::Vector3d &gpCoord,
                    const vector<double> &gpU,const vector<double> &gpV,
                    const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                    vector<double> &gpHist,const vector<double> &gpHistOld){
    
    switch (matetype){
        case MaterialType::NULLMATE:
            return;
        case MaterialType::CONSTPOISSONMATE:
            ConstPoissonMate(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        case MaterialType::NONLINEARPOISSONMATE:
            NonlinearPoissonMate(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        case MaterialType::CONSTDIFFUSIONMATE:
            ConstDiffusionMate(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        case MaterialType::NONLINEARDIFFUSIONMATE:
            NonlinearDiffusionMate(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        case MaterialType::CAHNHILLIARDMATE:
            CahnHilliardMate(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        case MaterialType::GENERALWAVEMATE:
            GenralWaveMate(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        case MaterialType::LINEARELASTICMATE:
            LinearElasticMate(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        case MaterialType::FINITESTRAINMATE:
            FiniteStrainMate(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        case MaterialType::STVENANTMATE:
            SaintVenantMate(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        case MaterialType::NEOHOOKEANMATE:
            NeoHookeanMate(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        case MaterialType::ELASTICTHERMALMEHCANICSMATE:
            ThermalElasticMate(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        case MaterialType::LINEARELASTICPHASEFIELDFRACTUREMATE:
            LinearElasticPhaseFieldFracMate(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        case MaterialType::MODIFYELASTICFRACTUREMATE:
            ModifyElasticFracMate(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        case MaterialType::USER1MATE:
            UserMate1(nDim,t,dt,_MaterialBlockList[iblock-1]._MateParams,gpCoord,gpU,gpV,gpGradU,gpGradV,gpHist,gpHistOld);
            break;
        default:
            Msg_MaterialSystem_UnsupportUmat();
            Msg_AsFem_Exit();
            break;
    }
}