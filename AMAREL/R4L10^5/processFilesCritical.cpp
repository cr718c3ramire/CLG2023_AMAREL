#include <stdio.h>
#include "auxCLG2023.h" 
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
 
/* CREATION DATE: MARCH 24, 2023
 * 
 * PURPOSE: MONTE CARLO BY READING FILES CONTAINING harmonic means averages in the forms of vectors. 
 *  
 * NOTES:
 * 
 * 0. THIS WAS A PROTOTYPE FOR A CLASS/HEADER FILE WHICH KEPT FAILING TO INITALIZE A VECTOR<VECTOR<DOUBLE>> 
 *    DATA MEMBER, SO I DID THIS PROCEDURAL VERSION
 * 
 * 1. I run this on the cluster since I don't want to download tons of files from cluster to my local PC. 
 
 * 
 * 
 * */ 

int main(int argc, char *argv[])
{
    ifstream File;
    vector<int>reportFileVector(7);
    long double number;
    

    const int pow2=27;
    int NUMBER_SUCCESSFUL_TASKS;//<---WE WILL CHECK IF THIS MATCHES NUMBERofTASKSinJOBARRAY, for safety
    int NUMBERofTASKSinJOBARRAY=10; 
    
    vector<vector<double>> HARMONIC_MEANS(NUMBERofTASKSinJOBARRAY, vector<double>(pow2));
     

   /*====================================================================================
    *  The following, given wildCard1, wildCard2:
    * 
      1. opens each SummVEC reportfile by using wildCard1 as a base name
          for each file which in a FOR loop appends the TASK_ID.
           * 
      2.  Once the TASK_ID is appended,the fileName is completely determined
           and we open the SummVEC file, check if successful run. 
      
      3. If successful, we use File2Name which depends on wildCard2+TASK_ID, and
          open the corresponding hmeansAvge file which has averaged the hmean of each 
          of NT trials run by this TASK_ID.  
    *==================================================================================== */   
   string numberTasksStr=to_string(NUMBERofTASKSinJOBARRAY);
   
   string parameters="R4L10^4_N8115_2^27_NT1_";        
   string wildCard1="SummVECR4L10^4_N8115_2^27_NT1_27182677_";
   string wildCard2="hmeanAvgeR4L10^4_N8115_2^27_NT1_27182677_";
   
   string file1Name=wildCard1;
   string file2Name=wildCard2;
   string taskIDstr="1";
   
   parameters=parameters+"on_"+numberTasksStr+"_tasks";
   
   /*====================================================================================
    THE ACTUAL FOR LOOP that loads data, WE DON'T HAVE TO INPUT ANYTHIING HERE. 
    *==================================================================================== */ 
   
   for(int taskNum=0; taskNum<NUMBERofTASKSinJOBARRAY; taskNum++)
   {
        
       taskIDstr=to_string(taskNum+1);//<--shift by 1 since AMAREL is 1-based index, not 0-based
       file1Name=wildCard1+taskIDstr;
       file2Name=wildCard2+taskIDstr;
       
       File.open(file1Name);
       for(int elem=0; elem< 7; elem++)
        {   
            File>>number;
            reportFileVector[elem]=number;
            }
       File.close();     
       
       if(reportFileVector[6]>0){cout<<"Opened a file of NT successful runs! \n"<<endl;}
       
       
       File.open(file2Name);
       for(int elem=0; elem< pow2; elem++)
        {   
            File>>number;
            HARMONIC_MEANS[taskNum][elem]=number;
            }
        
        File.close();    
            
        cout<<"We have so far at task number "<< taskNum<< ", HARMONIC_MEANS["<<taskNum<<"]= \n"<<endl;  
        DisplayVector(HARMONIC_MEANS[taskNum]);
        cout<<"\n";
       
       }
       
   /*====================================================================================
    THE AVERAGING DATA LOOP, also gives us a vector returning averages 
    *==================================================================================== */
    
    long double runningMean=0.0;
    
    vector<double> averageAcrossFiles(pow2);

    for(int hMeanNumber=0; hMeanNumber<pow2; hMeanNumber++)
    {
        runningMean=0.0;
        for(int taskNum=0; taskNum<NUMBERofTASKSinJOBARRAY; taskNum++)
        {
            number=HARMONIC_MEANS[taskNum][hMeanNumber];
            runningMean+=number;
            }
        runningMean/=NUMBERofTASKSinJOBARRAY;
        averageAcrossFiles[hMeanNumber]=runningMean;
    
            
        }
        
     cout<<"Obtained averageAcrossFiles= \n"<<endl;   
     DisplayVector(averageAcrossFiles);
     
     
   /*====================================================================================
        SAVE THE AVERAGE DATA ARRAY! 
    *==================================================================================== */
    
 
    
    ofstream avgeVecFile;
    string fileName="AvgeFileHM"+parameters; //<--eliminated date in FileName, could cause problms
    avgeVecFile.open(fileName);
    PythonFILEPrintVector(averageAcrossFiles, avgeVecFile);
    avgeVecFile.close();

}
    
