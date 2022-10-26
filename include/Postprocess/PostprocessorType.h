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
//+++ Date   : 2022.09.22
//+++ Purpose: defines postprocessor type in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

enum class PostprocessorType{
    NULLPPS,
    //*********************************
    // for nodal type pps
    //*********************************
    NODALVALUE,
    NODALSCALARMATERIALVALUE,
    NODALVECTORMATERIALVALUE,
    NODALRANK2MATERIALVALUE,
    NODALRANK4MATERIALVALUE,
    //*********************************
    // for side type pps
    //*********************************
    AREA,
    SIDEAVERAGEVALUE,
    SIDEAVERAGESCALARMATERIALVALUE,
    SIDEAVERAGEVECTORMATERIALVALUE,
    SIDEAVERAGERANK2MATERIALVALUE,
    SIDEAVERAGERANK4MATERIALVALUE,
    //
    SIDEINTEGRATEVALUE,
    SIDEINTEGRATESCALARMATERIALVALUE,
    SIDEINTEGRATEVECTORMATERIALVALUE,
    SIDEINTEGRATERANK2MATERIALVALUE,
    SIDEINTEGRATERANK4MATERIALVALUE,
    //
    SIDEMINVALUE,
    SIDEMINSCALARMATERIALVALUE,
    SIDEMINVECTORMATERIALVALUE,
    SIDEMINRANK2MATERIALVALUE,
    SIDEMINRANK4MATERIALVALUE,
    //
    SIDEMAXVALUE,
    SIDEMAXSCALARMATERIALVALUE,
    SIDEMAXVECTORMATERIALVALUE,
    SIDEMAXRANK2MATERIALVALUE,
    SIDEMAXRANK4MATERIALVALUE,
    //*********************************
    // for volume pps
    //*********************************
    VOLUME,
    //
    VOLUMEAVERAGEVALUE,
    VOLUMEAVERAGESCALARMATERIALVALUE,
    VOLUMEAVERAGEVECTORMATERIALVALUE,
    VOLUMEAVERAGERANK2MATERIALVALUE,
    VOLUMEAVERAGERANK4MATERIALVALUE,
    //
    VOLUMEINTEGRATEVALUE,
    VOLUMEINTEGRATESCALARMATERIALVALUE,
    VOLUMEINTEGRATEVECTORMATERIALVALUE,
    VOLUMEINTEGRATERANK2MATERIALVALUE,
    VOLUMEINTEGRATERANK4MATERIALVALUE,
    //
    MAXVALUE,
    MAXSCALARMATERIALVALUE,
    MAXVECTORMATERIALVALUE,
    MAXRANK2MATERIALVALUE,
    MAXRANK4MATERIALVALUE,
    //
    MINVALUE,
    MINSCALARMATERIALVALUE,
    MINVECTORMATERIALVALUE,
    MINRANK2MATERIALVALUE,
    MINRANK4MATERIALVALUE,
    //*********************************
    // for user-defined pps
    //*********************************
    USER1NODALPPS,
    USER2NODALPPS,
    USER3NODALPPS,
    USER4NODALPPS,
    USER5NODALPPS,
    //
    USER1SIDEINTEGRALPPS,
    USER2SIDEINTEGRALPPS,
    USER3SIDEINTEGRALPPS,
    USER4SIDEINTEGRALPPS,
    USER5SIDEINTEGRALPPS,
    //
    USER1VOLUMEINTEGRALPPS,
    USER2VOLUMEINTEGRALPPS,
    USER3VOLUMEINTEGRALPPS,
    USER4VOLUMEINTEGRALPPS,
    USER5VOLUMEINTEGRALPPS
};