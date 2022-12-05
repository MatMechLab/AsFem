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
//+++ Purpose: implement postprocess system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Utils/MessagePrinter.h"

#include "MathUtils/Vector3d.h"
#include "MathUtils/Rank2Tensor.h"
#include "MathUtils/Rank4Tensor.h"
#include "MathUtils/Vector.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "ProjectionSystem/ProjectionSystem.h"
#include "MateSystem/MateSystem.h"
#include "SolutionSystem/SolutionSystem.h"
#include "FE/FE.h"


#include "Postprocess/PostprocessorType.h"
#include "Postprocess/PostprocessorBlock.h"

/**
 * For different postprocessors
 */
#include "Postprocess/NodalValuePostprocessor.h"
#include "Postprocess/NodalScalarMatePostprocessor.h"
#include "Postprocess/NodalVectorMatePostprocessor.h"
#include "Postprocess/NodalRank2MatePostprocessor.h"
#include "Postprocess/NodalRank4MatePostprocessor.h"

// for side integral postprocessors
#include "Postprocess/AreaPostprocessor.h"
#include "Postprocess/SideIntegralValuePostprocessor.h"
#include "Postprocess/SideIntegralScalarMatePostprocessor.h"
#include "Postprocess/SideIntegralVectorMatePostprocessor.h"
#include "Postprocess/SideIntegralRank2MatePostprocessor.h"
#include "Postprocess/SideIntegralRank4MatePostprocessor.h"
// for user-X side integral postprocessors
#include "Postprocess/User1SideIntegralPostprocessor.h"

// for volume integral postprocessors
#include "Postprocess/VolumePostprocessor.h"
#include "Postprocess/VolumeIntegralValuePostprocessor.h"
#include "Postprocess/VolumeIntegralScalarMatePostprocessor.h"
#include "Postprocess/VolumeIntegralVectorMatePostprocessor.h"
#include "Postprocess/VolumeIntegralRank2MatePostprocessor.h"
#include "Postprocess/VolumeIntegralRank4MatePostprocessor.h"
#include "Postprocess/User1VolumeIntegralPostprocessor.h"

/**
 * This class implement the general postprocess (i.e., nodal value, side integration, volume integration) in AsFem
 */
class Postprocessor:public NodalValuePostprocessor,
                    public NodalScalarMatePostprocessor,
                    public NodalVectorMatePostprocessor,
                    public NodalRank2MatePostprocessor,
                    public NodalRank4MatePostprocessor,
                    // for side integral pps
                    public AreaPostprocessor,
                    public SideIntegralValuePostprocessor,
                    public SideIntegralScalarMatePostprocessor,
                    public SideIntegralVectorMatePostprocessor,
                    public SideIntegralRank2MatePostprocessor,
                    public SideIntegralRank4MatePostprocessor,
                    public User1SideIntegralPostprocessor,
                    // for volume integral pps
                    public VolumePostprocessor,
                    public VolumeIntegralValuePostprocessor,
                    public VolumeIntegralScalarMatePostprocessor,
                    public VolumeIntegralVectorMatePostprocessor,
                    public VolumeIntegralRank2MatePostprocessor,
                    public VolumeIntegralRank4MatePostprocessor,
                    public User1VolumeIntegralPostprocessor{
public:
    /**
     * constructor
     */
    Postprocessor();

    /**
     * initialize the postprocess class
     */
    void init();

    /**
     * set up the input filename
     * @param filename the string name of the input file
    */
   inline void setInputFileName(const string &filename){m_inputfilename=filename;}

    /**
     * add block to the list
     */
    void addPPSBlock2List(const PostprocessorBlock &t_block);

    /**
     * get the block number of pps
     */
    inline int getPPSBlocksNum()const{return m_pps_blocksnum;}
    /**
     * get the i-th pps block
     * @param i the block's index, starts from 1
     */
    inline PostprocessorBlock getIthPPSBlock(const int &i)const{
        if(i<1||i>m_pps_blocksnum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of pps blocks range");
            MessagePrinter::exitAsFem();
        }
        return m_pps_blocklist[i-1];
    }
    /**
     * get the i-th pps block's name
     * @param i the block's index, starts from 1
     */
    inline string getIthPPSBlockName(const int &i)const{
        if(i<1||i>m_pps_blocksnum){
            MessagePrinter::printErrorTxt("i="+to_string(i)+" is out of pps blocks range");
            MessagePrinter::exitAsFem();
        }
        return m_pps_blocklist[i-1].m_block_name;
    }
    /**
     * get the i-th pps block's name
     * @param i the block's index, starts from 1
     */
    inline string getCSVFileName()const{
        return m_csv_filename;
    }

    /**
     * execute the postprocess for different functions
     * @param t_mesh the mesh class
     * @param t_dofhandler the dofhandler class
     * @param t_fe the fe space
     * @param t_matesystem the material system class
     * @param t_projsystem the projection system class
     * @param t_solution the solution system class
     */
    void executePostprocess(const Mesh &t_mesh,
                            const DofHandler &t_dofhandler,
                            FE &t_fe,
                            MateSystem &t_matesystem,
                            ProjectionSystem &t_projsystem,
                            SolutionSystem &t_solution);

    /**
     * save all the postprocess results into csv file
     * @param time the current simulation time
     */
    void savePPSResults2CSVFile(const double &time);
    /**
     * write out the header info for csv file
     * @param inputfilename the string name of the input file
     */
    void prepareCSVFileHeader();

    /**
     * get the active status of the postprocess system
     */
    inline bool hasPostprocess()const{
        if(m_pps_blocksnum) return true;
        return false;
    }

    /**
     * get the output interval for postprocess output
    */
    inline int getInterval()const{
        return m_output_interval;
    }

    /**
     * release the allocated memory
     */
    void releaseMemory();

    /**
     * print out the postprocessor information
     */
    void printInfo()const;

private:
    /**
     * execute the nodal type postprocess
     * @param pps_type the type of postprocess
     * @param dofid the local id of the specific dof
     * @param t_parameters the json parameters taken from input file
     * @param t_dofhandler the dof handler class
     * @param t_soln the solution system
     * @param t_projsystem the projection system
     */
    double executeNodalPostprocess(const PostprocessorType &pps_type,
                                   const int &dofid,
                                   const nlohmann::json &t_parameters,
                                   const DofHandler &t_dofhandler,
                                   SolutionSystem &t_soln,
                                   ProjectionSystem &t_projsystem);
    /**
     * execute the side integral type postprocess
     * @param pps_type the type of postprocess
     * @param dofid the local id of the specific dof
     * @param sidenames the sides name vector
     * @param t_parameters the json parameters taken from input file
     * @param t_mesh the mesh class
     * @param t_dofhandler the dof handler class
     * @param t_fe the fe space class
     * @param t_soln the solution system
     * @param t_projsystem the projection system
     */
    double executeSideIntegralPostprocess(const PostprocessorType &pps_type,
                                          const int &dofid,
                                          const vector<string> &sidenames,
                                          const nlohmann::json &t_parameters,
                                          const Mesh &t_mesh,
                                          const DofHandler &t_dofhandler,
                                          FE &t_fe,
                                          SolutionSystem &t_soln,
                                          ProjectionSystem &t_projsystem);
    
    /**
     * execute the side integral postprocess libs
     * @param pps_type the type of postprocess
     * @param dofid the global id of the specific dof
     * @param nodeid the global id of the specific node
     * @param t_parameters the json parameters taken from input file
     * @param t_elmtinfo the local element info structure
     * @param t_shp the shape function class
     * @param t_soln the solution system
     * @param t_projsystem the projection system
     */
    double runSideIntegralPostprocessLibs(const PostprocessorType &pps_type,
                                          const int &dofid,
                                          const int &nodeid,
                                          const nlohmann::json &t_parameters,
                                          const LocalElmtInfo &t_elmtinfo,
                                          const LocalShapeFun &t_shp,
                                          SolutionSystem &t_soln,
                                          ProjectionSystem &t_projsystem);

    /**
     * execute the volume integral type postprocess
     * @param pps_type the type of postprocess
     * @param dofid the local id of the specific dof
     * @param domainnames the domain name vector
     * @param t_parameters the json parameters taken from input file
     * @param t_mesh the mesh class
     * @param t_dofhandler the dof handler class
     * @param t_fe the fe space class
     * @param t_soln the solution system
     * @param t_projsystem the projection system
     */
    double executeVolumeIntegralPostprocess(const PostprocessorType &pps_type,
                                            const int &dofid,
                                            const vector<string> &domainnames,
                                            const nlohmann::json &t_parameters,
                                            const Mesh &t_mesh,
                                            const DofHandler &t_dofhandler,
                                            FE &t_fe,
                                            SolutionSystem &t_soln,
                                            ProjectionSystem &t_projsystem);
    /**
     * execute the volume integral postprocess libs
     * @param pps_type the type of postprocess
     * @param dofid the global id of the specific dof
     * @param nodeid the global id of the specific node
     * @param t_parameters the json parameters taken from input file
     * @param t_elmtinfo the element info structure
     * @param t_shp the shape function class
     * @param t_soln the solution system
     * @param t_projsystem the projection system
     */
    double runVolumeIntegralPostprocessLibs(const PostprocessorType &pps_type,
                                            const int &dofid,
                                            const int &nodeid,
                                            const nlohmann::json &t_parameters,
                                            const LocalElmtInfo &t_elmtinfo,
                                            const LocalShapeFun &t_shp,
                                            SolutionSystem &t_soln,
                                            ProjectionSystem &t_projsystem);

private:
    string m_inputfilename;/** the string name of the input file*/
    string m_csv_filename;/**< the string name of the csv file */
    int m_output_interval;/**< the output interval */
    vector<string> m_pps_namelist;/**< the name vector for the final postprocessed variables */
    vector<double> m_pps_values;/**< the value vector for the final postprocessed variables */

    //***********************
    vector<PostprocessorBlock> m_pps_blocklist;/**< for the postprocess block defined in input file */
    int m_pps_blocksnum;/**< number of pps blocks */

    LocalElmtInfo m_local_elmtinfo;/**< for the local element info data structure */
    LocalShapeFun m_local_shp;/**< for the local shape function */
    Nodes m_nodes0,m_nodes;/**< the nodal coordinates of current element */
    Rank2Tensor m_xs;
    Vector3d m_normal;/**< the normal vector */

private:
    PetscMPIInt m_rank;/**< the local rank */
    PetscMPIInt m_size;/**< the number of total cpus */

};