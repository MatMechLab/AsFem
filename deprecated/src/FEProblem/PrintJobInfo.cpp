#include "FEProblem/FEProblem.h"

void FEProblem::PrintJobInfo(){
    cout<<"*** FE problem job information:                            ***"<<endl;
    printf("***    openmp threads=%4d                                 ***\n",Eigen::nbThreads());
    if(JobType==FEJobType::STATIC){
        cout<<"***    job type= Static                                    ***"<<endl;
    }
    else{
        if(feCtrlInfo._TimeSteppingMethod==TimeSteppingType::BACKWARDEULER){
            cout<<"***    job type= transient, method=backward euler(BE)      ***"<<endl;
        }
        else if(feCtrlInfo._TimeSteppingMethod==TimeSteppingType::CRANKNICOLSON){
            cout<<"***    job type= transient, method=cranknicolson(CN)       ***"<<endl;
        }
        else if(feCtrlInfo._TimeSteppingMethod==TimeSteppingType::BDF2){
            cout<<"***    job type= transient, method=BDF2                    ***"<<endl;
        }
        else if(feCtrlInfo._TimeSteppingMethod==TimeSteppingType::BDF3){
            cout<<"***    job type= transient, method=BDF3                    ***"<<endl;
        }
        else if(feCtrlInfo._TimeSteppingMethod==TimeSteppingType::BDF4){
            cout<<"***    job type= transient, method=BDF4                    ***"<<endl;
        }
        else if(feCtrlInfo._TimeSteppingMethod==TimeSteppingType::BDF5){
            cout<<"***    job type= transient, method=BDF5                    ***"<<endl;
        }
        printf("***    total time=%14.6e, dt0=%14.6e       ***\n",feCtrlInfo._EndTime,feCtrlInfo._dt0);
        if(feCtrlInfo._IsAdaptive){
            printf("***    adaptive=true,growcoef=%5.2f,cutcoef=%5.2f,iters=%2d ***\n",
               timeStepping.GetGrowthFactor(),timeStepping.GetCutbackFactor(),timeStepping.GetOptimIters());
       }
       else{
           cout<<"***    adaptive=false, using constant dt for stepping      ***"<<endl;
       }
       
    }
    if(feCtrlInfo._IsProjOn){
        cout<<"***    projection=true                                     ***"<<endl;
    }
    else{
        cout<<"***    projection=false                                    ***"<<endl;
    }
    cout<<"***--------------------------------------------------------***"<<endl;

}