#ifndef ASFEM_SOLUTIONSYSTEM_H
#define ASFEM_SOLUTIONSYSTEM_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cassert>


#include "Eigen/Eigen"
#include "Eigen/Core"

#include "DofHandler/DofHandler.h"

using namespace std;


class DofHandler;

class SolutionSystem
{
public:
    SolutionSystem();
    void SetDofsNum(int dofs) {_nDofs=dofs;}
    void SetNodesNum(int nodes) {_nNodes=nodes;}
    void SetElmtsNum(int nelmts) {_nBulkElmts=nelmts;}
    void SetProjsNum(int nproj) {_nProjsPerNode=nproj;}
    void SetHistsNum(int nhists) {_nHistValuePerGPoint=nhists;}
    void SetGPointsNumPerBulkElmt(int ngpoints) {_nGPointPerBulkElmt=ngpoints;}

    void SetProjNameFromStrVec(vector<string> &name){
        if(name.size()<1){
            cout<<"*** Error: projection name list is too short or empty !!!  ***"<<endl;
            cout<<"***        name=v1 v2... should be given !!!               ***"<<endl;
            Msg_AsFem_Exit();
        }
        _ProjNameList=name;_nProjsPerNode=int(_ProjNameList.size());
        _HasProjNameList=true;
    }

    inline int GetHistNumPerGPoint() const{
        return _nHistValuePerGPoint;
    }

    inline int GetProjNumPerNode() const{
        return _nProjsPerNode;
    }
    inline string GetIthProjName(const int &i) const{
        return _ProjNameList[i-1];
    }
    inline vector<string> GetProjNameList()const{
        return _ProjNameList;
    }

    void Init();

    void UpdateSolution(const int &currentstep);
    void UpdateHistValues(){_HistOld=_Hist;}

    inline double GetIthUNorm(const int &i) const {
        if(_UNorms.size()<1) return 0.0;
        if(i>(int)_UNorms.size()||i<1) assert("i is out of Unorms range in SolutionSystem.h !!!");
        return _UNorms[i-1];
    }

    inline double GetIthNodeJthDofValue(DofHandler &dofHandler,const int &i,const int &j) const{
        return _U(dofHandler.GetIthNodeJthDofFlag(i,j)-1);
    }
    inline double GetIthNodeJthOldDofValue(DofHandler &dofHandler,const int &i,const int &j) const{
        return _Uold(dofHandler.GetIthNodeJthDofFlag(i,j)-1);
    }
    inline double GetIthNodeJthOlderDofValue(DofHandler &dofHandler,const int &i,const int &j) const{
        return _Uolder(dofHandler.GetIthNodeJthDofFlag(i,j)-1);
    }
    //*****************************************************
    inline double GetIthNodeJthDofVeloValue(DofHandler &dofHandler,const int &i,const int &j) const{
        return _V(dofHandler.GetIthNodeJthDofFlag(i,j)-1);
    }
    inline double GetIthNodeJthOldDofVeloValue(DofHandler &dofHandler,const int &i,const int &j) const{
        return _Vold(dofHandler.GetIthNodeJthDofFlag(i,j)-1);
    }
    inline double GetIthNodeJthOlderDofVeloValue(DofHandler &dofHandler,const int &i,const int &j) const{
        return _Volder(dofHandler.GetIthNodeJthDofFlag(i,j)-1);
    }
    //***************************************************
    inline double GetIthNodeJthDofAcceValue(DofHandler &dofHandler,const int &i,const int &j) const{
        return _A(dofHandler.GetIthNodeJthDofFlag(i,j)-1);
    }
    inline double GetIthNodeJthOldDofAcceValue(DofHandler &dofHandler,const int &i,const int &j) const{
        return _Aold(dofHandler.GetIthNodeJthDofFlag(i,j)-1);
    }
    inline double GetIthNodeJthOlderDofAcceValue(DofHandler &dofHandler,const int &i,const int &j) const{
        return _Aolder(dofHandler.GetIthNodeJthDofFlag(i,j)-1);
    }
    //***************************************************
    // inline double GetIthElmtJthHistValue(const int &e,const int &j) const{
    //     return _Hist(e-1,(j-1)*_nGPointPerBulkElmt);
    // }
    // inline double& GetIthElmtJthHistValue(const int &e,const int &j){
    //     return _Hist((e-1)*_nHistValuePerElmt+j-1);
    // }
    // inline double GetIthElmtJthHistOldValue(const int &e,const int &j) const{
    //     return _HistOld((e-1)*_nHistValuePerElmt+j-1);
    // }
    // inline double& GetIthElmtJthHistOldValue(const int &e,const int &j){
    //     return _HistOld((e-1)*_nHistValuePerElmt+j-1);
    // }
    // inline void GetIthElmtHistValue(const int &e,vector<double> &localhist) const {
    //     if((int)localhist.size()>=_nHistValuePerElmt){
    //         for(int i=0;i<_nHistValuePerElmt;++i){
    //             localhist[i]=_Hist((e-1)*_nHistValuePerElmt+i);
    //         }
    //     }
    //     else{
    //         localhist.clear();
    //         for(int i=0;i<_nHistValuePerElmt;++i){
    //             localhist.push_back(_Hist((e-1)*_nHistValuePerElmt+i));
    //         }
    //     }
    // }
    // inline void GetIthElmtHistOldValue(const int &e,vector<double> &localhist) const {
    //     if((int)localhist.size()>=_nHistValuePerElmt){
    //         for(int i=0;i<_nHistValuePerElmt;++i){
    //             localhist[i]=_HistOld((e-1)*_nHistValuePerElmt+i);
    //         }
    //     }
    //     else{
    //         localhist.clear();
    //         for(int i=0;i<_nHistValuePerElmt;++i){
    //             localhist.push_back(_HistOld((e-1)*_nHistValuePerElmt+i));
    //         }
    //     }
    // }
    //***************************************************
    inline double GetIthNodeJthProjValue(const int &i,const int &j) const{
        return _ProjValue((i-1)*_nProjsPerNode+j-1);
    }
    inline double& GetIthNodeJthProjValue(const int &i,const int &j){
        return _ProjValue((i-1)*_nProjsPerNode+j-1);
    }

public:
    Eigen::VectorXd _U,_Uold,_Uolder; // for solution and historial solutions
    Eigen::VectorXd _V,_Vold,_Volder; // for velocity
    Eigen::VectorXd _A,_Aold,_Aolder; // for acceleration
    Eigen::VectorXd _Unew,_Utemp;// this is the one used for NR iterations

    Eigen::MatrixXd _ProjValue; // projected value of nodal point
    Eigen::MatrixXd _Hist,_HistOld; // history value of gauss point(elemental variable,not nodal one!!!)

private:
    vector<double> _UNorms;

private:
    int _nDofs,_nNodes,_nBulkElmts;
    int _nDofsPerNode;
    int _nProjsPerNode;
    int _nHistValuePerGPoint,_nGPointPerBulkElmt;

private:
    vector<string> _ProjNameList;
    bool _HasProjNameList=false;
};

#endif // ASFEM_SOLUTIONSYSTEM_H