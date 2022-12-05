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
//+++ Date   : 2022.08.12
//+++ Purpose: save results to different file format
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "OutputSystem/OutputSystem.h"

void OutputSystem::saveResults2File(const int &t_step,
                                    const Mesh &t_mesh,
                                    const DofHandler &t_dofHandler,
                                    SolutionSystem &t_solution,
                                    ProjectionSystem &t_projection){
    if(t_step<0){
        // negative for static result output
        m_outputfile_name=m_inputfile_name.substr(0,m_inputfile_name.size()-5);// remove '.json'
        m_outputfile_name+=".vtu";
    }
    else{
        // for nonzero value, the output file should be step-dependent
        std::ostringstream ss;
        ss<<std::setfill('0')<<std::setw(8)<<t_step;
        m_outputfile_name=m_inputfile_name.substr(0,m_inputfile_name.size()-5);// remove ".i" extension name
        m_outputfile_name = m_outputfile_name+"-"+ss.str()+ ".vtu";
    }
    
    if(m_fileformat==ResultFileFormat::VTU){
        VTUWriter::saveResults(m_outputfile_name,t_mesh,t_dofHandler,t_solution,t_projection);
    }
}
//***********************************************
//*** for pvd file
//***********************************************
void OutputSystem::savePVDHead(){
    MPI_Comm_rank(PETSC_COMM_WORLD, &m_rank);
    if(m_rank == 0){
        m_pvdfile_name=m_inputfile_name.substr(0,m_inputfile_name.size()-5)+".pvd";// remove ".json" extension name
        std::ofstream out;
        out.open(m_pvdfile_name,std::ios::out);
        if (!out.is_open()){
            MessagePrinter::printErrorTxt("can\'t create a new pvd file(="+m_pvdfile_name+")!, please make sure you have the write permission");
            MessagePrinter::exitAsFem();
        }
        out<<"<?xml version=\"1.0\"?>\n";
        out<<"<VTKFile type=\"Collection\" version=\"0.1\"\n"
             "         byte_order=\"LittleEndian\"\n"
             "         compressor=\"vtkZLibDataCompressor\">\n";
        out<<"<Collection>"<<endl;
        out.close();
    }
}
void OutputSystem::savePVDEnd(){
    MPI_Comm_rank(PETSC_COMM_WORLD, &m_rank);
    if(m_rank==0){
        std::ofstream out;
        out.open(m_pvdfile_name,std::ios::app|std::ios::out);
        if (!out.is_open()){
            MessagePrinter::printErrorTxt("can\'t open pvd file(="+m_pvdfile_name+")!, please make sure you have the write permission");
            MessagePrinter::exitAsFem();
        }
        out<<"</Collection>\n";
        out<<"</VTKFile>\n";
        out.close();
    }
}
void OutputSystem::savePVDResults(const double &current_time){
    MPI_Comm_rank(PETSC_COMM_WORLD, &m_rank);
    if(m_rank==0){
        std::ifstream in;
        string line;
        vector<string> lines;
        in.open(m_pvdfile_name,std::ios::in);
        while(!in.eof()){
            getline(in,line);
            lines.push_back(line);
        }
        in.close();

        // remove the last 3 lines
        lines.pop_back();
        lines.pop_back();
        lines.pop_back();

        std::ofstream out;
        char val[20];
        out.open(m_pvdfile_name,std::ios::out);
        if (!out.is_open()){
            MessagePrinter::printErrorTxt("can\'t open pvd file(="+m_pvdfile_name+")! please make sure you have the write permission");
            MessagePrinter::exitAsFem();
        }
        for(const auto &line:lines){
            out<<line<<endl;
        }

        sprintf(val,"%14.6e",current_time);

        out<<"<DataSet timestep=\""<<string(val)<<"\" "
           <<"group=\"\" part=\"0\" "
           <<"file=\""<<m_outputfile_name<<"\"/>"<<endl;
        out<<"</Collection>\n";
        out<<"</VTKFile>\n";
        out.close();
    }
}