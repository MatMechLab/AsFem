#include "MaterialSystem/MaterialSystem.h"

MaterialSystem::MaterialSystem()
{
    _nMaterialBlocks=0;
    _MaterialBlockList.clear();
    _nMateValues=50;
    _nRank2MateValues=6;
    _nRank4MateValues=3;
    _MateValues.clear();
    _Rank2TensorMateValues.clear();
    _Rank4TensorMateValues.clear();
}