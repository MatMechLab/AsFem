//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.05.12
//+++ Purpose: run the bulk element libs for the residual and jacobian
//+++          calculation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/BulkElmtSystem.h"

void BulkElmtSystem::runBulkElmtLibs(const FECalcType &t_calctype,const double (&ctan)[3],
                                     const int &t_subelmtid,
                                     const MaterialsContainer &t_materialscontainer_old,
                                     const MaterialsContainer &t_materialscontainer,
                                     const LocalElmtInfo &t_elmtinfo,
                                     const LocalElmtSolution &t_elmtsoln,
                                     const LocalShapeFun &t_shp,
                                     MatrixXd &K,
                                     VectorXd &R){
    switch (m_elmtblock_list[t_subelmtid-1].m_elmttype)
    {
    case ElmtType::LAPLACEELMT:
        LaplaceElement::computeAll(t_calctype,t_elmtinfo,ctan,t_elmtsoln,t_shp,t_materialscontainer_old,t_materialscontainer,K,R);
        break;
    case ElmtType::SCALARBODYSOURCEELMT:
        ScalarBodySourceElement::computeAll(t_calctype,t_elmtinfo,ctan,t_elmtsoln,t_shp,t_materialscontainer_old,t_materialscontainer,K,R);
        break;
    case ElmtType::POISSONELMT:
        PoissonElement::computeAll(t_calctype,t_elmtinfo,ctan,t_elmtsoln,t_shp,t_materialscontainer_old,t_materialscontainer,K,R);
        break;
    case ElmtType::DIFFUSIONELMT:
        DiffusionElement::computeAll(t_calctype,t_elmtinfo,ctan,t_elmtsoln,t_shp,t_materialscontainer_old,t_materialscontainer,K,R);
        break;
    case ElmtType::ALLENCAHNELMT:
        AllenCahnElement::computeAll(t_calctype,t_elmtinfo,ctan,t_elmtsoln,t_shp,t_materialscontainer_old,t_materialscontainer,K,R);
        break;
    case ElmtType::MECHANICSELMT:
        MechanicsElement::computeAll(t_calctype,t_elmtinfo,ctan,t_elmtsoln,t_shp,t_materialscontainer_old,t_materialscontainer,K,R);
        break;
    case ElmtType::CAHNHILLIARDELMT:
        CahnHilliardElement::computeAll(t_calctype,t_elmtinfo,ctan,t_elmtsoln,t_shp,t_materialscontainer_old,t_materialscontainer,K,R);
        break;
    case ElmtType::KOBAYASHIELMT:
        KobayashiElement::computeAll(t_calctype,t_elmtinfo,ctan,t_elmtsoln,t_shp,t_materialscontainer_old,t_materialscontainer,K,R);
        break;
    case ElmtType::STRESSDIFFUSIONELMT:
        StressDiffusionElement::computeAll(t_calctype,t_elmtinfo,ctan,t_elmtsoln,t_shp,t_materialscontainer_old,t_materialscontainer,K,R);
        break;
    case ElmtType::MIEHEFRACTUREELMT:
        MieheFractureElement::computeAll(t_calctype,t_elmtinfo,ctan,t_elmtsoln,t_shp,t_materialscontainer_old,t_materialscontainer,K,R);
        break;
    case ElmtType::ALLENCAHNFRACTUREELMT:
        AllenCahnFractureElement::computeAll(t_calctype,t_elmtinfo,ctan,t_elmtsoln,t_shp,t_materialscontainer_old,t_materialscontainer,K,R);
        break;
    case ElmtType::STRESSCAHNHILLIARDELMT:
        StressCahnHilliardElement::computeAll(t_calctype,t_elmtinfo,ctan,t_elmtsoln,t_shp,t_materialscontainer_old,t_materialscontainer,K,R);
        break;
    case ElmtType::DIFFUSIONACFRACTUREELMT:
        DiffusionACFractureElement::computeAll(t_calctype,t_elmtinfo,ctan,t_elmtsoln,t_shp,t_materialscontainer_old,t_materialscontainer,K,R);
        break;
    default:
        MessagePrinter::printErrorTxt("unsupported bulk element type in runElmtLibs, please check your code");
        MessagePrinter::exitAsFem();
        break;
    }
}