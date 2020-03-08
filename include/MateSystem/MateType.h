//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2019
//* https://github.com/walkandthinker/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_MATETYPE_H
#define ASFEM_MATETYPE_H

enum class MateType{
    NullMate,
    ConstWaveMate,
    ConstDiffusionMate,
    NonLinearDiffusionMate,
    TensorDiffusionMate,
    ConstPoissonMate,
    NonLinearPoissonMate,
    TensorPoissonMate,
    LinearElasticMate,
    SaintVenantMate,
    FiniteStrainMate,
    NeoHookeanMate,
    MooneyRivlinMate,
    CahnHilliardMate,
    TensorCahnHilliardMate,
    LinearThermalMechanicsMate,
    LinearElasticCahnHilliardMate,
    NeoHookeanCahnHilliardMate,
    LinearElasticPhaseFieldFracMate,
    MieheLinearElasticFracMate,
    CohesivePFFracMate,
    MieheNeoHookeanFracMate,
    BordenLinearElasticFracMate,
    NeoHookeanPhaseFieldFracMate,
    User1Mate,
    User2Mate,
    User3Mate,
    User4Mate,
    User5Mate,
    User6Mate,
    User7Mate,
    User8Mate,
    User9Mate,
    User10Mate,
    User11Mate,
    User12Mate,
    User13Mate,
    User14Mate,
    User15Mate,
    User16Mate,
    User17Mate,
    User18Mate,
    User19Mate,
    User20Mate
};


#endif // ASFEM_MATETYPE_H