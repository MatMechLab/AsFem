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
//+++ Date   : 2022.09.23
//+++ Purpose: implement postprocess system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocessor.h"

void Postprocessor::prepareCSVFileHeader(){
    m_csv_filename=m_inputfilename.substr(0,m_inputfilename.size()-5)+".csv";
    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);
    if(m_rank==0){
        std::ofstream out;
        out.open(m_csv_filename.c_str(),std::ios::out);
        if(!out.is_open()){
            MessagePrinter::printErrorTxt("can\'t open "+m_csv_filename+", please make sure you have the write permission");
            MessagePrinter::exitAsFem();
        }
        out<<"time";
        for(auto it:m_pps_namelist){
            out<<","<<it;
        }
        out<<endl;
        out.close();
    }
}
void Postprocessor::savePPSResults2CSVFile(const double &time){
    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);
    if(m_rank==0){
        std::ofstream out;
        out.open(m_csv_filename.c_str(),std::ios::app|std::ios::out);
        if(!out.is_open()){
            MessagePrinter::printErrorTxt("can\'t open "+m_csv_filename+", please make sure you have the write permission");
            MessagePrinter::exitAsFem();
        }
        out<<std::scientific<<std::setprecision(8);
        out<<time;
        for(auto it:m_pps_values){
            out<<","<<it;
        }
        out<<endl;
        out.close();
    }
}