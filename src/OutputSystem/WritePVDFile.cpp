//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.02.20
//+++ Purpose: write results to pvd file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "OutputSystem/OutputSystem.h"

void OutputSystem::WritePVDFileHeader(){
    MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);
    if(_rank == 0){
        _PVDFileName=_InputFileName.substr(0,_InputFileName.size()-2)+".pvd";// remove ".i" extension name
        ofstream out;
        out.open(_PVDFileName,ios::out);
        if (!out.is_open()){
            string str="can\'t create a new pvd file(="+_PVDFileName+")!, please make sure you have write permission";
            MessagePrinter::PrintErrorTxt(str);
            MessagePrinter::AsFem_Exit();
        }
        out<<"<?xml version=\"1.0\"?>\n";
        out<<"<VTKFile type=\"Collection\" version=\"0.1\"\n"
             "         byte_order=\"LittleEndian\"\n"
             "         compressor=\"vtkZLibDataCompressor\">\n";
        out<<"<Collection>\n";
        out.close();
    }
}
//**********************************************
void OutputSystem::WritePVDFileEnd(){
    MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);
    if(_rank==0){
        ofstream out;
        out.open(_PVDFileName,ios::app|ios::out);
        if (!out.is_open()){
            string str="can\'t open pvd file(="+_PVDFileName+")!, please make sure you have write permission";
            MessagePrinter::PrintErrorTxt(str);
            MessagePrinter::AsFem_Exit();
        }
        out<<"</Collection>\n";
        out<<"</VTKFile>\n";
        out.close();
    }
}
//**********************************************
void OutputSystem::WriteResultToPVDFile(const double &timestep,string resultfilename){
    MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);
    if(_rank==0){
        ofstream out;
        char val[20];
        out.open(_PVDFileName,ios::app|ios::out);
        if (!out.is_open()){
            string str="can\'t open pvd file(="+_PVDFileName+")!, please make sure you have write permission";
            MessagePrinter::PrintErrorTxt(str);
            MessagePrinter::AsFem_Exit();
        }
        sprintf(val,"%14.6e",timestep);
        out<<"<DataSet timestep=\""<<string(val)<<"\" "
           <<"group=\"\" part=\"0\" "
           <<"file=\""<<resultfilename<<"\"/>"<<endl;
    }
}