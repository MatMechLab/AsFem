#ifndef ASFEM_OUTPUTSYSTEM_H
#define ASFEM_OUTPUTSYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

#include "Utils/StringUtils.h"
#include "MessagePrint/MessagePrint.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "SolutionSystem/SolutionSystem.h"

#include "Eigen/Eigen"

using namespace std;

class OutputSystem
{
public:
    OutputSystem();
    inline string GetLogFileName() const{return _LogFileName;}
    inline string GetVTUFileName() const{return _VTUFileName;}
    void InitOutputStream();
    void SetInputFileName(string inputname) {_InputFileName=inputname;_IsInit=false;}

    void WriteResultToVTU(Mesh &mesh,DofHandler &dofHandler,const Eigen::VectorXd &U);
    void WriteResultToVTU(Mesh &mesh,DofHandler &dofHandler,const Eigen::VectorXd &U,const int &nproj,const vector<string> &namelist,const Eigen::MatrixXd &Proj);
    
    // for time step output
    void WriteResultToVTU(const int &step,Mesh &mesh,DofHandler &dofHandler,const Eigen::VectorXd &U);
    void WriteResultToVTU(const int &step,Mesh &mesh,DofHandler &dofHandler,const Eigen::VectorXd &U,const int &nproj,const vector<string> &namelist,const Eigen::MatrixXd &Proj);

private:
    bool _IsInit=false;
    bool _IsLogOn=false;
    string _InputFileName;
    string _OutputFilePrefix;
    string _LogFileName;
    ofstream _VTUFile,_LogFile;
    string _VTUFileName;
};

#endif // ASFEM_OUTPUTSYSTEM_H