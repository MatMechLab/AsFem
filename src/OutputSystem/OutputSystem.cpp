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
//+++ Date   : 2020.07.12
//+++ Purpose: define the output system for AsFem, where all the 
//+++          results should be written out to the result file
//+++          by this class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "OutputSystem/OutputSystem.h"

OutputSystem::OutputSystem(){
    m_inputfile_name.clear();
    m_outputfile_name.clear();
    m_pvdfile_name.clear();
    m_fileformat=ResultFileFormat::VTU;
    m_intervals=1;
}

void OutputSystem::printInfo()const{
    MessagePrinter::printNormalTxt("Output information summary:");
    if(m_fileformat==ResultFileFormat::VTU){
        MessagePrinter::printNormalTxt("  output file format = vtu");
    }
    else if(m_fileformat==ResultFileFormat::VTK){
        MessagePrinter::printNormalTxt("  output file format = vtk");
    }
    else if(m_fileformat==ResultFileFormat::CSV){
        MessagePrinter::printNormalTxt("  output file format = csv");
    }
    MessagePrinter::printNormalTxt("  output interval = "+to_string(m_intervals));
    MessagePrinter::printStars();
}
