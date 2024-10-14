/* harmonic_mean_ladder.h
 * 
 * UPDATE 03/03/2023: I CREATED THIS AROUND AUGUST 2021 IT SEEMS! AND I HAVE BEEN CARRYING IT AROUND SINCE
 *ON THIS DATE I COPY-PASTED THIS HEADER FILE INTO MY CLG2023 FILE PROTOTYPE 
 * */
#ifndef HARMONIC_MEAN_LADDER_H
#define HARMONIC_MEAN_LADDER_H



/* ****************************  HARMONIC MEANS VECTOR IN BLOCKS OF LENGTH 2^j  ****************************
 * 
 * This Harmonic_Mean2k takes in a whole rho_a vector of size ULS=2^{k+1}-1, for k an integer. VERY IMPORTANT
 * THAT ULS=2^{k+1}-1 for this to work. For instance, if k=3, you would have the 2^4-1=15 values of rho_a with indices
 * divided into k+1 containers:
 *                
 *                      || 0 || 1 2 || 3 4 5 6 || 7 8 9 10 11 12 13 14 ||
 *   where the jth (j=0,1,...,k) container has 2^j elements. 

       Note that ULS-1=14=2^4-2, not 2^4-1. This will give k+1 harmonic means (the first one, 1/rho_a(0) is trivial   */ 

             
//template<typename rhoClass>
vector<double> Harmonic_Means2k(vector<double> & rho_a, int k)
{
    int size_means=k+1;
    vector<double> harmonic_means(size_means);
    int lower=1;
    int upper=2;
    int to_add=2; 
    double running_sum_inverses=0;
    double tempry=0;
    
    
    harmonic_means[0]=rho_a[0]; //notice rho_a size ULS>> k+1
    
    for(int j=1; j<size_means; ++j)
    {
        
        running_sum_inverses=0;
        //cout<<"note that lower=" <<lower<<" while upper= "<<upper<<endl;
        for(int idx=lower; idx<upper+1;++idx)
        {
            tempry=1/(rho_a[idx]);
            running_sum_inverses+=tempry;
            }
            harmonic_means[j]=to_add/running_sum_inverses; //notice at j of outer loop, to_add=2^j. This faster than pow(2,j) I
// believe  
            //cout<<"got it for j=" <<j<<endl;
            lower=upper+1;
            to_add=2*to_add;
            upper+=to_add;
            
        }
        
        return harmonic_means;
    
    
    }

/*******************************************************************************************

 *             ONE HARMONIC MEAN FUNCTION: I.E., RETURNS A SCALAR, NOT A VECTOR
                         
HARMONIC MEANS SIMILAR/LIKE MY MATLAB VERSION: THIS ONE ONLY RETURNS ONE HARMONIC_MEAN ON AN INTERVAL [S1, S2]
 * THE PREVIOUS HARMONIC_MEANS2K RETURNS A VECTOR, BUT ALSO, IT COMPUTES A HARMONIC MEAN IN BLOCKS OF LENGTH 
 * 2^j IN A MORE EFFICIENT MANNER THAN THE BELOW FUNCTION COULD DO IF INVOKED TO MIMIC THE ABOVE. THE ABOVE JUST KEEPS
 * TRACK MORE EASILY OF 2^j POWERS AS j=0,1,....,k. 
 * 
 

*********************************************************************************************/

// This function takes in a rho_a vector and calculates a harmonic mean from component s1 to component s2

template<typename rhoClass>
double Harmonic_Mean_Rho(int S1, int S2, rhoClass& rho_a) //computes harmonic mean on [S1, S2] (inclusive)
{
    double harmonic_mean=0;
    int length=S2-S1+1;
    double sum_inverses=0;
    double tempry=0;
    
    for(int s=S1;s<S2; ++s)  
    {
        tempry=1/rho_a[s];
        sum_inverses+=tempry;
        }
        
        harmonic_mean=double(length)/sum_inverses;
    return harmonic_mean;
    }


/* The reason I suggest two types of harmonic means functions is because with one we get say a vector double of harmonic
means. But this is in blocks size 2^j, the last one being a mean length 2^k, which is 2^k/ 2^(k+1)-1\approx 1/2 of the
total run info. It might be easier to look at the run, say if from 3000 to 50,000 it is easily in the good region, we use that 
 harmonic mean from [3,000, 50,000], instead of looking for some multiple of 2^j to add up */
 
 /*      ***********************   ***********************  ***********************   ***********************
  *                         FUNCTION THAT AVERAGES OVER A MATRIX
  *      ***********************   ***********************  ***********************   ***********************
  * */
  
  template<typename anArray>
  vector<double> AveragesByColumn(anArray& infoMatrix)
  {
      
      int numberOfTrials=infoMatrix.size(); /*notice the way my program does trials, trial number NT's harmonic means are
                                          in the NT-th row of infoMatrix. infoMatrix.size() returns only number of columns  */
      
      int numberOfColumns=infoMatrix[0].size(); 
      vector<double> resultVector(numberOfColumns);  //the ith entry stores averages of ith column
                                    
      double running_sum=0;

     for(int column=0;column< numberOfColumns; ++column)
     {
         running_sum=0;
         for(int row=0;row<numberOfTrials;++row)
         {
             running_sum+=infoMatrix[row][column];
             }
             
             resultVector[column]=running_sum/double(numberOfTrials);
         
         }
          
          return resultVector;                           
      }
  


#endif
