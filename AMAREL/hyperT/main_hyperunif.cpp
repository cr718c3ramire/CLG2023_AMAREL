#include <stdio.h>
#include<time.h>
#include<fstream>
#include<iostream>
#include<random>
#include<algorithm>
#include<iterator>
#include<array>
#include<cmath>
#include<string>
#include <sstream>
#include <chrono>
 
using namespace std::chrono;
using namespace std;

#include "CLG2023.h" 
#include "hyperunif2023.h"
#include "auxCLG2023.h"




/******************************************************************************************************************** 
 * Created April 3, 2023 for hyperuniform CLG2023
 ********************************************************************************************************************* */


int main(int argc, char *argv[])
{
    /* ---------------------------------------------------------------------------
     * Some auxiliary variables 
     *--------------------------------------------------------------------------- */
    auto start = high_resolution_clock::now();
    int r,i;
     
    

    cout<<"Starting program Rule 4 CLG2023 Hyperuniform"<<endl;
  
    
    /*---------------------------------------------------------------------------
     * DATA PARAMETERS SPECIFIC TO THIS RUN
     * 
     * To give in idea, for R4, L=10^4 at N=8100-8200 does well with 
     * ULS=2^27-1 steps, so pow2=27
     * 
     * for L=10^3, pow2=19 for N=820 particles
     * 
     * 
     * ---------------------------------------------------------------------------*/
    long double steps_this_trial {0};
    const int pow10=4;
    const size_t base=10;
   
    const size_t L=pow(base, pow10);
    size_t N=80;
    int NT=10;
    int pow2=29;
    vector<double> hmean{0};
    size_t ULS=pow(2,pow2)-1;
    int upperLimitFails=10000;
    
    
    /*------------------------------------------------------------------------------------ 
     * HYPERUNIFORM ASSISTANCE VARIABLES AND CONTAINERS
     *------------------------------------------------------------------------------------  */
    long double rho=N/float(2*L);

    vector<int> lengths{1,3,5, 7, 9, 15, 25, 35, 55,75, 95, 125, 225, 351, 425, 625, 1025, 1125, 1325, 1525, 1775,  2025, 2499};
    int num_lengths=lengths.size();
    vector<long double> averageVP2(num_lengths);
    vector<vector<long double>> VP2Data(NT, vector<long double>(num_lengths));
    //vector<vector<long double>> VP2Data[NT][num_lengths]{0}; //<--This won't be returned by functionm averageVP2 above is who matters
 
    int NUMBER_ICs_tried{0};
   

    int numberFrozenSoFar=0;
    vector<int> lifeSpan(NT);

    
    int running_nact{0};
    
    
    /****************************************************************************************************************************
             JUST TO GET SOME FILENAMES AND PARAMETERS FOR THEM, THIS ONLY FOR SLURM JOB ARRAYS. NOTE:
              * 
              * 1. RANDOM SEED SET USING TIME OF RUN AND THEN OFFSETTING THE TIME BY THE SLURM_ARRAY_TASK_ID 
              *    SO THAT EVEN IF JOBS START AT THE SAME TIME (WHICH HAPPENS OFTEN!) WE GET DIFFERENT RANDOM
              *    NUMBERS FOR EACH SLURM_ARRAY_TASK_ID
              * 
              * 2. BECAUSE IT REQUESTS getenv(SLURM_ARRAY_TASK_ID), I DON'T THINK THIS SCRIPT WILL MAKE SENSE
                   ON A SINGLE CORE/NON JOB ARRAY, SINCE I UNDERSTAND SLURM JOB ARRAYS ASSIGN SLURM_ARRAY_TASK_ID
 ****************************************************************************************************************************/

    /*
     * ADDED MARCH 19, 2023: A REPORT VECTOR: This will contain basically a vector form of Summary/Report
     * 
     *         reportVec=[NUMBER_TASK_SLURM_ARRAY, NT, L, N, ULS, time_took, SUCCESSFUL?]
     *                                       (SUCCESSFUL=1 if all NT ran well)
     * 
     * THIS REPORT VECTOR ADDED AT BOTTOM ONCE EVERYTHING RAN WELL..
     * 
     * PENDING: AT SOME POINT MAYBE A LIFESPAN ARRAY, BUT IT IS COMPLICATED BECAUSE REPORT VECTOR HAS  
     * "AVERAGED DATA" WHEREAS TO BUILD A LIFESPAN HISTOGRAM TYPE CHART, WE NEED THE RESULT OF EACH
        INDIVIDUAL TRIAL...
     * 
     * SHOULD BE A WHOLE LIFESPAN ARRAY AND THE MCFILES.H HEADER SHOULD JUST PUT THEM ALL TOGETHER AND 
     * THAT'S IT!
     * 
     * */
     
    vector<double> summryVector(7);
     
  
    string function_name="R4hyperU_";
    
    
    string stringL = to_string(pow10); //<---note this expresses L=10^power, such as 10^6
    string stringN=to_string(N);
    
    string stringNT=to_string(NT);
    string parameters="L10^"+stringL+"_"+"N"+stringN+"_NT"+stringNT; 
    
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, sizeof(buffer),"%m-%d-%Y",timeinfo);
    std::string time_str(buffer);
    
    
    parameters=parameters+'_'+time_str;
    
    
    long double seedRandom=time(0);
    string baseFileName=function_name+parameters+"_";
    
    char const* tmp = getenv( "SLURM_ARRAY_JOB_ID" );
    if ( tmp == NULL ) {
    //  Big problem...
         } else {
          std::string s( tmp );
          baseFileName+=s;
         // cout<< "We have " << s<<endl;
         }

    baseFileName+="_";
    
    char const* tmp2 = getenv( "SLURM_ARRAY_TASK_ID" );
    if ( tmp2 == NULL ) {
    //  Big problem...
         } else {
          std::string s2( tmp2 );
          baseFileName+=s2;
          seedRandom+=stod(s2); //<--OFFSET HERE USING TASK_ID
          summryVector[0]=stod(s2);
         }

    mt19937 mt_rand(seedRandom); //<-- using offset TASK_ID, get a random number
    auto rand_number=mt_rand();
    
    
    summryVector[1]=NT;
    summryVector[2]=L;
    summryVector[3]=N;
    summryVector[4]=ULS; //<--for hyperunif this not too relevant
      
  
      
 /****************************************************************************************************************************/
    
    
    
    /* ---------------------------------------------------------------------------
     * INITIALIZE MEMBERS OF OUR CLASS, AND OBJECT OF CLASS LightLadder 
     *--------------------------------------------------------------------------- */
     
    NbrsTable extDefTable(L);
    LightLadder<L> a;
    
    /*  *******************************************************************************************************************       
     *                              START THE ACTION! 
     *             
     * *******************************************************************************************************************/
     
    
   while(numberFrozenSoFar<NT){
    
    
    cout<<"Running " <<function_name<< " at L="   <<L<<", N="<<N<< ", NT= "<< NT<<" at ULS="<<ULS<<"\n"<<endl;
    a.steps_taken=0;
    a.ZeroOut();
    running_nact=0;
    steps_this_trial=0; 
    a.Populate(N, mt_rand); //<---populate array here with N 1's
    a.R4WhoActiveIC(extDefTable); //Check who's active in IC
    //cout<< "Obtained eta= "<<'\n';
    //cout<<"Running a new IC"<<endl; 
    a.ICR4noH(ULS, extDefTable, mt_rand);
    running_nact=a.GetNact();
    steps_this_trial=a.steps_taken;
    
    NUMBER_ICs_tried=NUMBER_ICs_tried+1;
    if(running_nact==0){
         cout<< "DEATH at step= "<<steps_this_trial<<"\n"<<endl;
         lifeSpan[numberFrozenSoFar]=steps_this_trial;
         //cout<< "numberFrozenSoFar="<<numberFrozenSoFar<<'\n';
         
         // Get data of this frozen config!
         HU(a.eta, lengths, VP2Data, rho, numberFrozenSoFar, L, NT);
         
         numberFrozenSoFar++;  
    }
    else{
        cout<<"Got a survivor, we will try again"<<'\n'; //<---this program does NOT capture statistics of surviving trial....maybe we should?? 
        }
    cout<<"Running RULE 4 L=" <<L<<" N="<<N<< " NT= "<< NT<<" at ULS="<<ULS<<"\n"<<endl;
    cout<<"Number of ICs attempted= "<< NUMBER_ICs_tried<<" while number ICs frozen= "<< numberFrozenSoFar <<"\n"<<endl;
     
    cout<<"We will end program if we accumulate number of FROZENS="<<NT<<'\n';
	 
    }
    
    
    AvgeHU(VP2Data, averageVP2, NT, num_lengths);
    
    
    /*  *******************************************************************************************************************       
     *                              END THE ACTION! FINAL PROCESSING/SAVING/ETC NOW BELOW:
     * 
     * MARCH 19, 2023: USED C++ FRIENDLY PRINTING TO FILE THIS TIME, INSTEAD OF PYTHON FRIENDLY:

                         VectorToFileCPP() 
     *             
     * *******************************************************************************************************************/

       
        ////////// STOP TIMING HERE!
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(stop - start);
        int secondsDuration=duration.count();
        cout <<"Program time took (rounded to seconds): "<< secondsDuration << endl;
        
         
        
        
        /* ********************************************************************************************************
         * ADDED MARCH 19, 2023: ZEROTH FILE: REPORT VECTOR SUMMARY TO MAKE PROCESSING DATA EASIER WHEN RUNNING
                                  ON MULTIPLE CORES 
         * 
         * reportVec=[NUMBER_TASK_SLURM_ARRAY, NT, L, N, ULS, runtime, SUCCESSFUL?]
         * 
         *                                                              ^^^SUCCESSFUL=1 if all NT ran well
         ********************************************************************************************************* */
         
        summryVector[5]=secondsDuration;
        summryVector[6]=1; //<---if we got up to this line, we are successful..
    
    
     
         
       
    
        
        ofstream lifeSpanFile;
        string lifeSpanFileName="lSpan"+ baseFileName;
        lifeSpanFile.open(lifeSpanFileName);
        VectorToFileCPP(lifeSpan, lifeSpanFile);
        lifeSpanFile.close();
        
       
    
        /*       SECOND FILE: REPORT FILE        */
        ofstream reportFile;
        string reportFileName="Summ"+  baseFileName; 
        reportFile.open(reportFileName);
        VectorToFileCPP(summryVector, reportFile);
        
        reportFile<<" \n\n DATE: "<<time_str<<"\n\n"<<endl;
       // reportFile<<"Program time took (rounded to seconds): "<< duration.count() << endl;
        reportFile<<" Program Rule 4 for hyperuniform at L="<<L<<" , N="<<N<<" ,NT=" <<NT<<" for ULS="<<ULS<<endl;
        reportFile<<" with NUMBER_ICs_tried (to get NT trials): "<< NUMBER_ICs_tried <<'\n'<<endl;
        reportFile<<" We accumulated data from lifespans of total fails accumulated: "<< NT <<'\n'<<endl;
        
    
        reportFile<<" The average (across NT trials) is array averageVP2=";
        PythonFILEPrintVector(averageVP2, reportFile);
        
        reportFile <<"\n the lengths used for this hyperuniform data were: \n"<<endl;
        PythonFILEPrintVector(lengths, reportFile);
 
        reportFile<<'\n';
        reportFile<<"\n \n Program up to here took a total time (rounded to nearest second): "<<secondsDuration<<endl;
        cout<<"Finished writing to file, will close now! "<<endl;
        reportFile.close();
        
        
        /*    THIS ONE THE CRUCIAL TO AVERAGE ACROSS AMAREL CORES: 3rd FILE CONTAINS averageVP2 array  

 * 
 *  
                            */
        ofstream avgeVP2File;
        string avgeVP2FileName="avgeVP2"+baseFileName;//<--eliminated date in FileName, could cause problms
        avgeVP2File.open(avgeVP2FileName);
        VectorToFileCPP(averageVP2, avgeVP2File);
         
        avgeVP2File.close();
        
        
     
         
    
    
    
    
}
    

 