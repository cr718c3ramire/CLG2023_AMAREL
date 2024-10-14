

/* CREATION DATE: 03/03/2023
 * PURPOSE: TO TEST MY CLG2023.H LIBRARY WHICH I STARTED TO DEVELOP LATE FEBRUARY 2023
 * 
 * 
 * 
 * */


#include <stdio.h>
#include "CLG2023.h" 
#include "harmonic_mean_ladder.h"
#include<time.h>
#include<fstream>
#include<iostream>
#include<random>
#include<algorithm>
#include<iterator>
#include<array>
#include<cmath>
#include<string>
#include<sstream>
#include<chrono>

#include<iomanip>
#include<ctime>
#include<sstream>
using namespace std;
 
using namespace std::chrono;
 
 

int main(int argc, char *argv[])
{
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
   
    const size_t ell=pow(base, pow10);
    size_t N=8110;
    int NT=1;
    int pow2=5;
    vector<double> hmean{0};
    long double ULS=pow(2,pow2)-1;
    int upperLimitFails=10000;
 
    /* ---------------------------------------------------------------------------
     * Some auxiliary variables and containers
     *--------------------------------------------------------------------------- */
    auto start = high_resolution_clock::now();
    int r,i;
    
       
    
    vector<vector<double>> HARMONIC_MEANS(NT, vector<double>(pow2)); //element harmonic_means[trial][interval]
    int NUMBER_ICs_tried{0};
    int NUMBER_ICs_survived{0};
    vector<int> deathStats(upperLimitFails);
    char ch;
    
    
    size_t running_nact{0};
    
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
     
  
    string function_name="R4";
    
    
    string stringL = to_string(pow10); //<---note this expresses L=10^power, such as 10^6
    string stringN=to_string(N);
    string stringIVALS=to_string(pow2); //so ran to ULS=2^pow2-1
    string stringNT=to_string(NT);
    string parameters="L10^"+stringL+"_"+"N"+stringN+"_"+"2^"+stringIVALS+"_"+"NT"+stringNT; 
    
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
    summryVector[2]=ell;
    summryVector[3]=N;
    summryVector[4]=ULS;
      
  
      
 /****************************************************************************************************************************/
    
    
    
    /* ---------------------------------------------------------------------------
     * INITIALIZE MEMBERS OF OUR CLASS, AND OBJECT OF CLASS LightLadder 
     *--------------------------------------------------------------------------- */
     
    NbrsTable extDefTable(ell);
    LightLadder<ell> a;
	 
     
    
    cout<<"Starting program "<< function_name << ":" <<endl;
    
    /*  *******************************************************************************************************************       
     *                              START THE ACTION! 
     *             
     * *******************************************************************************************************************/
    
    
    while(NUMBER_ICs_survived<NT){
    
    
           
    cout<<"Running " <<function_name<< " at L="   <<ell<<", N="<<N<< ", NT= "<< NT<<" at ULS="<<ULS<<"\n"<<endl;
   
    running_nact=0;
     
    a.Populate(N, mt_rand); //<---populate array here with N 1's
    a.R4WhoActiveIC(extDefTable); //Check who's active in IC
    //cout<< "Obtained eta= "<<'\n';
     
    hmean=a.OneICR4H2(pow2, extDefTable, mt_rand);
    running_nact=a.GetNact();
 
    
    NUMBER_ICs_tried=NUMBER_ICs_tried+1;
    if(running_nact==0){
         //DO NOTHING, run another IC

         cout<< "DEATH at step= "<<steps_this_trial<<"\n"<<endl;
         ;
    }
    else{
        HARMONIC_MEANS[NUMBER_ICs_survived]=hmean; 
        NUMBER_ICs_survived=NUMBER_ICs_survived+1;
        }
    a.ZeroOut(); //<---Flush it    
    cout<<"Running L=" <<ell<<" N="<<N<< " NT= "<< NT<<" at ULS="<<ULS<<"\n"<<endl;
    cout<<"Number of ICs attempted= "<< NUMBER_ICs_tried<<"\n"<<endl;
    cout<<"Number of ICs survived= "<< NUMBER_ICs_survived<<"\n"<<endl;
	
    }
    

    
    /*       HARMONIC_MEANS_AVERAGES[j] = harmonic average jth intval       */
    vector<double> HARMONIC_MEANS_AVERAGES(pow2);
    HARMONIC_MEANS_AVERAGES=AveragesByColumn(HARMONIC_MEANS);
    
    //Get sample standard deviation of HMEANS across NT
    
    vector<double> STDvector(pow2);
    double running_sum_for_std=0;
    double temporary=0;
    double tempHmean=0;
    
    for(int hmeanNum=0;hmeanNum<pow2;hmeanNum++)
    {
        tempHmean=HARMONIC_MEANS_AVERAGES[hmeanNum];
        running_sum_for_std=0;
        temporary=0;
        for(int j=0; j<NT; j++)
        {
            temporary=tempHmean- HARMONIC_MEANS[j][hmeanNum];
            running_sum_for_std+=temporary*temporary;
            }
            
        running_sum_for_std=running_sum_for_std/double(NT-1);  
        STDvector[hmeanNum]=sqrt(running_sum_for_std);
        
        }
    
    
    
    
    
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
        
        time_t rawtime;
        struct tm * timeinfo;
        char buffer[80];

        time (&rawtime);
        timeinfo = localtime(&rawtime);

        strftime(buffer,sizeof(buffer),"%m-%d-%Y",timeinfo);
        std::string time_str(buffer);
        
        
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
         
        ofstream summVecFile;
        string summVecFileName="SummVEC"+baseFileName; //<--eliminated date in FileName, could cause problms
        summVecFile.open(summVecFileName);
        VectorToFileCPP(summryVector, summVecFile);
        summVecFile.close();
        
        
        
    
        /*       SECOND FILE: REPORT FILE CONTAINS HARMONIC MEAN OVER LARGE TAIL
                 and human-legible summary of run. Note above summVEC is a vector version of this report
                 that is easy for the computer to read and further process.   
                          
                           * */
        ofstream reportFile;
        string reportFileName="Summ"+baseFileName+ '_' +time_str;
        reportFile.open(reportFileName);
       // reportFile<<"Program time took (rounded to seconds): "<< duration.count() << endl;
        reportFile<<" DATE: "<<time_str<<"\n\n"<<endl;
         
         
        
        reportFile<<"Program " <<function_name<< " at L="   <<ell<<", N="<<N<< ", NT= "<< NT<<" at ULS="<<ULS<<"\n"<<endl;
        reportFile<<" with NUMBER_ICs_tried (to get NT trials): "<< NUMBER_ICs_tried <<'\n'<<endl;
        
        reportFile<<'\n'<< "The standard deviation array from each harmonic mean is: "<< '\n'<< '\n';
        PythonFILEPrintVector(STDvector,reportFile); //<--not too useful
    
        reportFile<<"\n \n Program up to here took a total time (rounded to nearest second): "<<secondsDuration<<endl;
        reportFile<<"\n The vector form of the summary here is: \n"<<endl;
        PythonFILEPrintVector(summryVector,reportFile);//<---also printed in first file above
        cout<<"Finished writing to file, will close now! "<<endl;
        
        reportFile.close();
        
        /*    THIS ONE TO PLOT: 3rd FILE CONTAINS HARMONIC_MEANS_AVERAGES[j]= avge harmonic mean across trials on jth interval length 2^j   */
        ofstream harmonicAveragesFile;
        string hmeanFileName="hmeanAvge"+baseFileName;//<--eliminated date in FileName, could cause problms
        harmonicAveragesFile.open(hmeanFileName);
        VectorToFileCPP(HARMONIC_MEANS_AVERAGES, harmonicAveragesFile);
         
        harmonicAveragesFile.close();
        
        
         
    
    
    
    
}
