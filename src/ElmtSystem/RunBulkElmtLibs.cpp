//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.11.29
//+++ Purpose: we list all of our elements here for different models
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/BulkElmtSystem.h"

void BulkElmtSystem::RunBulkElmtLibs(const FECalcType &calctype,const ElmtType &elmtytype,
                         const double (&ctan)[2],
                         const LocalElmtInfo &elmtinfo,
                         const LocalElmtSolution &soln,
                         const LocalShapeFun &shp,
                         const Materials &Mate,const Materials &MateOld,
                         ScalarMateType &gpProj,
                         MatrixXd &localK,VectorXd &localR){
    switch (elmtytype){
        case ElmtType::LAPLACEELMT:
            break;
        case ElmtType::POISSONELMT:
            PoissonElmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::TIMEDERIVELMT:
            break;
        case ElmtType::DIFFUSIONELMT:
            DiffusionElmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::CAHNHILLIARDELMT:
            CahnHilliardElmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::MECHANICSELMT:
            MechanicsElmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::MIEHEFRACELMT:
            MieheFractureElmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER1ELMT:
            User1Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER2ELMT:
            User2Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER3ELMT:
            User3Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER4ELMT:
            User4Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER5ELMT:
            User5Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER6ELMT:
            User6Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER7ELMT:
            User7Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER8ELMT:
            User8Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER9ELMT:
            User9Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER10ELMT:
            User10Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER11ELMT:
            User11Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER12ELMT:
            User12Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER13ELMT:
            User13Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER14ELMT:
            User14Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER15ELMT:
            User15Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER16ELMT:
            User16Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER17ELMT:
            User17Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER18ELMT:
            User18Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER19ELMT:
            User19Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        case ElmtType::USER20ELMT:
            User20Elmt::ComputeAll(calctype,elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj,localK,localR);
            break;
        default:
            MessagePrinter::PrintErrorTxt("unsupported element type in ElmtSystem, please check your code or your input file");
            MessagePrinter::AsFem_Exit();
            break;
    }
}
