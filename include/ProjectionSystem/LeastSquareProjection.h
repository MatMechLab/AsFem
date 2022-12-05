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
//+++ Date   : 2022.07.22
//+++ Purpose: Implement least square projection from integration
//+++          points to nodal points
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "ProjectionSystem/ProjectionBase.h"

/**
 * This class implements the least square projection from integration points to nodal points
 */
class LeastSquareProjection:public ProjectionBase{
public:
    LeastSquareProjection();

protected:
    /**
     * the global projection action
     * @param t_mesh the mesh class
     * @param t_data the projection data structure
     */
    virtual void globalProjectionAction(const Mesh &t_mesh,
                                        ProjectionData &t_data) override;
    /**
     * the local projection action
     * @param nodesnum the nodes number of current element
     * @param t_elconn the local element connectivity
     * @param detjac the jacobian determinte of current qpoint
     * @param t_shp the shapefunction class
     * @param m_mate the material container from material system
     * @param t_data the projection data structure
     */
    virtual void localProjectionAction(const int &nodesnum,
                                       const vector<int> &t_elconn,
                                       const double &detjac,
                                       const ShapeFun &t_shp,
                                       const MaterialsContainer &t_mate,
                                       ProjectionData &t_data) override;

private:
    /**
     * the local projection action for scalar material
     * @param nodesnum the nodes number of current element
     * @param elconn the elemental connectivity, start from 1
     * @param detjac the jacobian determinte of current qpoint
     * @param t_shp the shapefunction class
     * @param t_mate the material container from material system
     * @param nproj the number of scalar material to be projected
     * @param t_scalar_name the scalar material name list
     * @param t_scalar_projvec the scalar projection vector
     */
    void projectLocalScalarMate2Global(const int &nodesnum,
                                       const vector<int> &elconn,
                                       const double &detjac,
                                       const ShapeFun &t_shp,
                                       const MaterialsContainer &t_mate,
                                       const int &nproj,
                                       const vector<string> &t_scalar_name,
                                       Vector &t_scalar_projvec);
    
    /**
     * the local projection action for vector material
     * @param nodesnum the nodes number of current element
     * @param elconn the elemental connectivity, start from 1
     * @param detjac the jacobian determinte of current qpoint
     * @param t_shp the shapefunction class
     * @param t_mate the material container from material system
     * @param nproj the number of vector material to be projected
     * @param t_vector_name the vector material name list
     * @param t_vector_projvec the vector projection vector
     */
    void projectLocalVectorMate2Global(const int &nodesnum,
                                       const vector<int> &elconn,
                                       const double &detjac,
                                       const ShapeFun &t_shp,
                                       const MaterialsContainer &t_mate,
                                       const int &nproj,
                                       const vector<string> &t_vector_name,
                                       Vector &t_vector_projvec);
    /**
     * the local projection action for rank-2 material
     * @param nodesnum the nodes number of current element
     * @param elconn the elemental connectivity, start from 1
     * @param detjac the jacobian determinte of current qpoint
     * @param t_shp the shapefunction class
     * @param t_mate the material container from material system
     * @param nproj the number of rank-2 material to be projected
     * @param t_rank2_name the rank-2 material name list
     * @param t_rank2_projvec the rank-2 projection vector
     */
    void projectLocalRank2Mate2Global(const int &nodesnum,
                                      const vector<int> &elconn,
                                      const double &detjac,
                                      const ShapeFun &t_shp,
                                      const MaterialsContainer &t_mate,
                                      const int &nproj,
                                      const vector<string> &t_rank2_name,
                                      Vector &t_rank2_projvec);

    /**
     * the local projection action for rank-4 material
     * @param nodesnum the nodes number of current element
     * @param elconn the elemental connectivity, start from 1
     * @param detjac the jacobian determinte of current qpoint
     * @param t_shp the shapefunction class
     * @param t_mate the material container from material system
     * @param nproj the number of rank-4 material to be projected
     * @param t_rank4_name the rank-4 material name list
     * @param t_rank4_projvec the rank-4 projection vector
     */
    void projectLocalRank4Mate2Global(const int &nodesnum,
                                      const vector<int> &elconn,
                                      const double &detjac,
                                      const ShapeFun &t_shp,
                                      const MaterialsContainer &t_mate,
                                      const int &nproj,
                                      const vector<string> &t_rank4_name,
                                      Vector &t_rank4_projvec);

};
