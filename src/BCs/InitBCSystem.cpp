//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "BCs/BCSystem.h"

void BCSystem::InitBCSystem(Mesh &mesh){
    _xi=0.0;_eta=0.0;_JxW=0.0;_nDim=0.0;_nNodesPerBCElmt=0;
    _nDim=0;
    _nNodesPerBCElmt=mesh.GetNodesNumPerBulkElmt();
    _elNodes.InitNodes(mesh.GetNodesNumPerBulkElmt());
}