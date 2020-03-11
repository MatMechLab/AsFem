//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_ELMTTYPE_H
#define ASFEM_ELMTTYPE_H

enum class ElmtType{
    NullElmt,
    PoissonElmt,
    DiffusionElmt,
    WaveElmt,
    CahnHilliardElmt,
    TensorCahnHilliardElmt,
    MechCahnHilliardElmt,
    MechanicsElmt,
    ThermalElmt,
    ThermalMechanicsElmt,
    PhaseFieldFracElmt,
    MieheFractureElmt,
    CohesivePFFracElmt,
    BordenFractureElmt,
    DendriteElmt,
    User1Elmt,
    User2Elmt,
    User3Elmt,
    User4Elmt,
    User5Elmt,
    User6Elmt,
    User7Elmt,
    User8Elmt,
    User9Elmt,
    User10Elmt,
    User11Elmt,
    User12Elmt,
    User13Elmt,
    User14Elmt,
    User15Elmt,
    User16Elmt,
    User17Elmt,
    User18Elmt,
    User19Elmt,
    User20Elmt
};

#endif //ASFEM_ELMTTYPE_H