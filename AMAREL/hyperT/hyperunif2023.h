#ifndef HYPERUNIF2023_H
#define HYPERUNIF2023_H
 
#include<iostream>
 
#include<algorithm>
#include<iterator>
#include "auxCLG2023.h"
using namespace std;

/*
 * CREATED  MARCH 26, 2023 BASED ON main_hyperuniform.cpp within Rule4CriticalD project from August 2022 
 * 
 * THIS TO WORK WITH CLG2023 WHERE ARRAY IS 1D!
 * 
 * 
 * 
 * */
 
  
     
  /*The result should be in vector<long double> averageVP2(num_lengths), BUT MODIFIED BY REFERENCE!;
   * 
   * THUS THE AIM IS TO:
   * GOAL: MODIFY averageVP2 ARRAY BY REFERENCE
   * 
   * TO ACHIEVE GOAL IT NEEDS TO BE FED ALSO: 
   * 1. ETA=FINAL FROZEN CONFIG
   * 2. KNOW NT: IF FROM JOBARRAY, THEN      
   * 
   * I'm assuming lengths can be fed into this function as 
   * 
   * vector<int> lengths{1,2, 3, 4, 5, 7, 8,  9, 13, 16, 25, 32, 35, 55, 64, 75,
                 95, 115, 128, 225, 256, 425, 512, 625, 1024, 1125, 1525, 1775,  2048, 2499,
                 3000, 3335, 4096, 5025, 7125, 8192, 12125, 16384, 21125, 32768, 51125, 65536,
                 83335, 97765, 115005, 131072, 175075, 225125, 333333};
   * 
   * 
   * */
  
  
 template<typename etaClass, typename vectorClass1, typename vectorClass2>
 void HU(const etaClass& eta,  vectorClass1& lengths, vectorClass2& VP2Data,
 long double rho, int numberFrozenSoFar, size_t L, int NT)
 {
              
     int num_lengths=lengths.size();
      

     /*SOME TEMPRY VARIABLES TO ASSIST*/
    
    int L_over_l=0;
    int ell=0;
    int N_i_1=0;
    int N_i_2=0;
    int N_i=0;
    int colSubVolLower=0;
    int colSubVolUpper=0;
 
    long double runningSumVP2=0;
    long double temp1=0;
    long double temp2=0;
    long double temp3=0;
    
    /*NOTE THAT 1-D ARRAY USED OFTEN FOR MY ETA ARRAYS. THUS NEED TO USE 2DTO1D CONVERTING TYPE FUNCTIONS */
    
    size_t row=0;
    size_t col=0;
    
    for(int ellNumber=0; ellNumber<num_lengths; ellNumber++)
           {
          /*---------------------------------------------
           * EXTRACT SOME N_k's
           *--------------------------------------------- */
               ell=lengths[ellNumber];
              
               L_over_l=L/ell;
              
               colSubVolLower=0;
               colSubVolUpper=ell;
               
               
               runningSumVP2=0;
                
               for(int subVolNum=0; subVolNum<L_over_l; subVolNum++)
               {
                   N_i_1=0;
                   N_i_2=0;
                   for(int col=colSubVolLower; col<colSubVolUpper; col++)
                   {
                       N_i_1+=eta[Index1D(0,col)];
                       N_i_2+=eta[Index1D(1,col)];
                       }
    
                   N_i=N_i_1+N_i_2;    
                   
                   temp3=(N_i/float(2*ell)-rho);
                   runningSumVP2+=temp3*temp3;    
                   
                   colSubVolLower=colSubVolUpper;
                   colSubVolUpper+=ell;
                   
                   
                   }
              
               VP2Data[numberFrozenSoFar][ellNumber]=runningSumVP2/double(L_over_l);
                
               }
               
     
     }
 
 
 
template<typename VP2DataClass, typename vectorClass1>
 void AvgeHU(VP2DataClass& VP2Data, vectorClass1& averageVP2, int NT, int num_lengths)
 { 
     long double temp=0;
     for(int ellNumber=0; ellNumber<num_lengths; ellNumber++)
           {
               averageVP2[ellNumber]=0;
               temp=0;
                

               for(int trialNumber=0;trialNumber<NT; trialNumber++)
               {
                   temp+=VP2Data[trialNumber][ellNumber];
                   }

               temp=temp/NT;
               averageVP2[ellNumber]=temp;
                
               }
     
 }
#endif