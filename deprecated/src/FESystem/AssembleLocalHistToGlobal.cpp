#include "FESystem/FESystem.h"

void FESystem::AssembleLocalHistToGlobal(const int &elmtid,const int &gpInd,const vector<double> &gpHist,Eigen::MatrixXd &Hist){
    for(int i=1;i<=_nHist;++i){
        Hist.coeffRef(elmtid-1,(gpInd-1)*_nHist+i-1)=gpHist[i-1];
    }
}