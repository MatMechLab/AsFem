//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_OUTPUTFILETYPE
#define ASFEM_OUTPUTFILETYPE

enum class OutputFileType{
    NullType,
    VTK,
    VTU,
    CSV,
    TXT,
    VTP
};

#endif