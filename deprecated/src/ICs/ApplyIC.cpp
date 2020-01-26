#include "ICs/ICSystem.h"

void ICSystem::ApplyIC(Mesh &mesh,DofHandler &dofHandler,Eigen::VectorXd &U){
    int iblock,DofIndex;
    BCBlock icblock;
    for(iblock=1;iblock<=GetICBlockNum();++iblock){
        DofIndex=dofHandler.GetDofIndexViaName(_ICBlockList[iblock-1]._DofName);
        switch(_ICBlockList[iblock-1]._ICTypeID){
            case ICType::CONST:
                ApplyConstIC(_ICBlockList[iblock-1]._value,DofIndex,mesh,dofHandler,U);
                break;
            case ICType::RAND:
                ApplyRandIC(_ICBlockList[iblock-1]._value,DofIndex,mesh,dofHandler,U);
                break;
            case ICType::CIRCLE:
                ApplyCircleIC(_ICBlockList[iblock-1]._value,DofIndex,mesh,dofHandler,U);
                break;
            default:
                cout<<"*** Error: unsupported IC type!!!                     ***"<<endl;
                Msg_AsFem_Exit();
                break;
        }
    }
}