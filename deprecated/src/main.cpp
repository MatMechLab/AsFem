#include <iostream>


//* For AsFem's own header file
#include "Welcome.h"
#include "Eigen/Eigen"
#include "FEProblem/FEProblem.h"

// #include "Utils/RankTwoTensor.h"
// #include "Utils/RankFourTensor.h"

int main(int args,char *argv[])
{
    Eigen::initParallel();
    const double AsFem_Version=0.2;
    const int AsFem_Year=2019;
    const int AsFem_Mon=11;
    const int AsFem_Day=30;
    Welcome(AsFem_Version,AsFem_Year,AsFem_Mon,AsFem_Day);

    FEProblem feProblem(args,argv);
    
    feProblem.RunFEProblem();




    return 0;
}
