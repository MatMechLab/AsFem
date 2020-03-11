//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::RunElmtLib(const int &isw,const ElmtType &iuel,const int &nDim,const int &nNodes,
                const double &t,const double &dt,const double (&ctan)[2],
                const Vector3d &gpCoord,
                const vector<double> &gpU,const vector<double> &gpV,
                const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
                const ShapeFun &shp,
                const vector<double> &ScalarMaterials,
                const vector<Vector3d> &VectorMaterials,
                const vector<RankTwoTensor> &Rank2Materials,
                const vector<RankFourTensor> &Rank4Materials,
                vector<double> &gpHist,const vector<double> &gpHistOld,vector<double> &gpProj,
                MatrixXd &localK,VectorXd &localR){
    switch (iuel){
    case ElmtType::NullElmt:
        return;
        break;
    case ElmtType::PoissonElmt:
        Poisson(isw,nDim,nNodes,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                shp,
                ScalarMaterials,VectorMaterials,Rank2Materials,Rank4Materials,
                gpHist,gpHistOld,gpProj,localK,localR);
        break;
    case ElmtType::DiffusionElmt:
        Diffusion(isw,nDim,nNodes,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                shp,
                ScalarMaterials,VectorMaterials,Rank2Materials,Rank4Materials,
                gpHist,gpHistOld,gpProj,localK,localR);
        break;
    case ElmtType::CahnHilliardElmt:
        CahnHilliard(isw,nDim,nNodes,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                shp,
                ScalarMaterials,VectorMaterials,Rank2Materials,Rank4Materials,
                gpHist,gpHistOld,gpProj,localK,localR);
        break;
    case ElmtType::MechCahnHilliardElmt:
        MechanicalCahnHilliard(isw,nDim,nNodes,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                shp,
                ScalarMaterials,VectorMaterials,Rank2Materials,Rank4Materials,
                gpHist,gpHistOld,gpProj,localK,localR);
        break;
    case ElmtType::TensorCahnHilliardElmt:
        TensorCahnHilliard(isw,nDim,nNodes,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                shp,
                ScalarMaterials,VectorMaterials,Rank2Materials,Rank4Materials,
                gpHist,gpHistOld,gpProj,localK,localR);
        break;
    case ElmtType::MechanicsElmt:
        Mechanics(isw,nDim,nNodes,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                shp,
                ScalarMaterials,VectorMaterials,Rank2Materials,Rank4Materials,
                gpHist,gpHistOld,gpProj,localK,localR);
        break;
    case ElmtType::ThermalMechanicsElmt:
        ThermalMechanics(isw,nDim,nNodes,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                shp,
                ScalarMaterials,VectorMaterials,Rank2Materials,Rank4Materials,
                gpHist,gpHistOld,gpProj,localK,localR);
        break;
    case ElmtType::MieheFractureElmt:
        MieheFracture(isw,nDim,nNodes,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                shp,
                ScalarMaterials,VectorMaterials,Rank2Materials,Rank4Materials,
                gpHist,gpHistOld,gpProj,localK,localR);
        break;
    case ElmtType::CohesivePFFracElmt:
        CohesivePFFracture(isw,nDim,nNodes,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                shp,
                ScalarMaterials,VectorMaterials,Rank2Materials,Rank4Materials,
                gpHist,gpHistOld,gpProj,localK,localR);
        break;
    case ElmtType::BordenFractureElmt:
        BordenFracture(isw,nDim,nNodes,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                shp,
                ScalarMaterials,VectorMaterials,Rank2Materials,Rank4Materials,
                gpHist,gpHistOld,gpProj,localK,localR);
        break;
    case ElmtType::DendriteElmt:
        DendriteModel(isw,nDim,nNodes,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                shp,
                ScalarMaterials,VectorMaterials,Rank2Materials,Rank4Materials,
                gpHist,gpHistOld,gpProj,localK,localR);
        break;
    case ElmtType::ThermalElmt:
        Thermal(isw,nDim,nNodes,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                shp,
                ScalarMaterials,VectorMaterials,Rank2Materials,Rank4Materials,
                gpHist,gpHistOld,gpProj,localK,localR);
        break;
    case ElmtType::User1Elmt:
        UserElmt1(isw,nDim,nNodes,t,dt,ctan,gpCoord,gpU,gpV,gpGradU,gpGradV,
                shp,
                ScalarMaterials,VectorMaterials,Rank2Materials,Rank4Materials,
                gpHist,gpHistOld,gpProj,localK,localR);
        break;
    default:
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported element type in uel                      !!!   ***\n");
        Msg_AsFem_Exit();
        break;
    }
}