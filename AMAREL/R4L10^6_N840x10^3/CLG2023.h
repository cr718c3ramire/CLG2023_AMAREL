#ifndef CLG2023_H
#define CLG2023_H
#include "auxCLG2023.h"
#include<iostream>
#include<cstddef>
#include<random>
#include<algorithm>
#include<iterator>
using namespace std;

/* APPARENTLY, CAN ONLY IMPLEMENT TEMPLATE CLASSES WITHIN THE HEADER FILE, UNLESS SOME
 * OTHER THING IS DONE. I THINK [STEFANO] USES SOME "NAMESPACE" NOTATION SO HE CAN 
 * IMPLEMENT IN A SEPARATE .CPP FILE. I DIDN'T WANT TO DO THAT...*/

/****************************************************************************************
 * Feb 27, 2023: I am seriously thinking of just keeping the nbrs table as an object or
 * separate class but done without template size_t L so it isn't built over and over 
 * 
 * PLUS I JUST REMEMBERED THERE IS ALMOST NO POINT IN TEMPLATE<SIZE_T> UNLESS IT WILL REALLY 
 * SAVE METHODS/FUNCTIONS RELYING ON L WITH LOTS OF LOOPS. THIS TABLE NBRS ONLY BUILT ONCE,
 * SO NO NEED TO TEMPLATE CLASS IT. JUST GIVE IT A CONSTRUCTOR LENGTH L NbrTable nbrs(L);
 * 
 * Then a la [Joshi] just store a reference pointer to that nbrs object and dereference. 
 * 
 * Feb 27, 2023 log #2: Ended up storing reference to NbrsTable object. 
 *
 ***************************************************************************************** */
 
 
class NbrsTable
{

  public:
      size_t** elems;
      size_t tLength;
      NbrsTable(size_t ell):tLength{ell}, elems{new size_t*[2*ell]}
      {
          //cout<<"Creating nbrs table: \n" <<endl;
          for(size_t i=0; i<2*ell; i++){elems[i]=new size_t[3];}
          for(size_t i=2; i<2*ell-2; i++)
          {
              elems[i][0]=i-2;
              elems[i][1]=i+2;
              
              if(IsOdd(i)){elems[i][2]=i-1;}   
              else{elems[i][2]=i+1;}

              }
          
          elems[0][0]=2*tLength-2;
          elems[0][1]=2;
          elems[0][2]=1;
         
          elems[1][0]=2*tLength-1;
          elems[1][1]=3;
          elems[1][2]=0;
         
          elems[2*tLength-2][0]=2*tLength-4;
          elems[2*tLength-2][1]=0;
          elems[2*tLength-2][2]=2*tLength-1;
         
          elems[2*tLength-1][0]=2*tLength-3;
          elems[2*tLength-1][1]=1;
          elems[2*tLength-1][2]=2*tLength-2;    
              
 
      }
      ~NbrsTable()
      {
          for(size_t i=0; i<2*tLength; i++)
          {
              //cout<< "Destroying in NbrsTable Class for each site i = "<< i<<" the ptr:\n"<<endl;
              delete[] elems[i];
                   }
           //cout<< "Destroying in NbrsTable Class the most outer ptr:\n"<<endl;        
           delete[] elems;       
              
          } 
 
 
};
 
/****************************************************************************************
 * ****************************************************************************************
THE ACTUAL WORKHORSE CLASS. I IMPLEMENT SOME OF THE THINGS WITHIN IT, JUST THE 
 CONSTRUCTOR AND DESTRUCTOR AND SOME MINOR INLINE OPERATORS THAT ARE USEFUL
  * 
  * NOTE I USE SOME AUX FUNCTIONS FROM auxCLG2023.H
  * 
****************************************************************************************
****************************************************************************************/ 

 
template<size_t L>
class LightLadder
{
private: 
     const size_t length;
     size_t NACT; //<---GetNact() public method returns NAC
     size_t* active_sts; //<--1D index notation. Returns 1-D index of active_sts
     size_t* pos;//<--a la Speer
     int* empty; //<--a la Speer
      
     
     /**********************************************************************************
      * SOME DUMB LOCAL VARIABLES ONLY USEFUL SO I DON'T REDEFINE OVER AND OVER  
      * AT EACH OF POTENTIALLY 10^9 MOVES
      * 
      **********************************************************************************/
      
     
     int* status_updates;
     int number_locals_changed=0;
     
     size_t old_across_nbr=0;
     size_t old_left_nbr=0;
     size_t old_right_nbr=0;

     size_t local_sites_changed[6]{0};
     
     int eta_left=0;
     int eta_across=0;
     int eta_right=0;
   
    /**********************************************************************************/
public:
     int* eta;
     //const NbrsTable & nbrs; //<--AMarel compiler wants NbrsTable const &!
     size_t siteSelected;
     size_t steps_taken=0; 
      
     
 
     LightLadder():length{L}, eta{new int[2*L]}, 
     active_sts{new size_t[2*L]}, pos{new size_t[2*L]}, empty{new int[2*L]}, NACT{0}, 
     siteSelected{0}, status_updates{new int[2*L]}
     {
         //cout<< "Creating a LightLadder: \n"<<endl;
         for(size_t idx =0; idx<2*L; idx++)
             {
                 eta[idx]=0;
                 active_sts[idx]=0;
                 pos[idx]=0;
                 empty[idx]=0;
                 status_updates[idx]=0;
             }

     }
     
     ~LightLadder()
     {
         delete[] eta;
         delete[] active_sts;
         delete[] pos;
         delete[] empty;
         delete[] status_updates;
         
         }
     
     inline int& operator()(size_t x, size_t y) {return eta[Index1D(x,y)]; }
     inline const int& operator()(size_t x, size_t y) const { return eta[Index1D(x,y)]; }
     
     /*-------------------------------------------------------------------------
      * THE ACTUAL EVOLUTION RULES, IMPLEMENTED ELSEWHERE
      *----------------------------------------------------------------------*/
     void R2OneMove(const NbrsTable& nbrs);
     void R4OneMove(const NbrsTable& nbrs);
     vector<double> OneICR4H2(const int power2, const NbrsTable& nbrs,  mt19937& mt);
     
    
     /*-------------------------------------------------------------------------
      * SOME AUXILIARY  NOT SO TRIVIAL FUNCTIONS.
      * 
      *----------------------------------------------------------------------*/
     size_t GetLength() const {return length;}
     size_t GetNact() const {return NACT;}
     
     void R4WhoActiveIC(const NbrsTable& nbrs); //<--fill in active_sts using R4's defn of active or not
      
     void Populate(size_t N, mt19937& mt);//<--you'll be fed the mt by driver program
     void ZeroOut();
     
     /*-------------------------------------------------------------------------
      * SOME TRIVIAL FUNCTIONS
      *----------------------------------------------------------------------*/
      
      void PrintEta() const;
  
   
   
};

/***********************************************************************************/
/***********************************************************************************/
/*================================================================================== 
 * IMPLEMENTATIONS OF CLASS LightLadder<L> 
 *================================================================================== */
 /***********************************************************************************/
 /***********************************************************************************/
 
template<size_t L> 
void LightLadder<L>::PrintEta() const
{
    for(int r=0;r<2; r++)
            {
                for(size_t i=0; i<L; i++)
                {
                  cout<< eta[Index1D(r,i)]<< '\t';
                }
        
                 cout << '\n';
            }
            cout << '\n';
    
}

template<size_t L> 
void LightLadder<L>::ZeroOut() 
{
    //cout<<"We Zero out: \n" <<endl;
    for(size_t idx =0; idx<2*L; idx++)
    {
        eta[idx]=0;
        pos[idx]=0;
        active_sts[idx]=0;
        empty[idx]=0;
        }
    
    }

/*
 * The following Populate() uses RandomSample1() in auxCLG2023.h, mersenne twister mt fed by driver 
 * program main()
 * 
 * */

template<size_t L> 
void LightLadder<L>::Populate(size_t N, mt19937& mt) 
 {
     //cout<<"Calling .Populate() method: \n" <<endl;
     vector<size_t> indices_to_shuffle(2*L); 
     for (size_t i=0; i<2*L; i++){ indices_to_shuffle[i]=i;}
 
     shuffle (indices_to_shuffle.begin(), indices_to_shuffle.end(), mt);

     vector<size_t> indices_chosen(N);
     for(size_t i=0;i<N; i++){indices_chosen[i]=indices_to_shuffle[i];}
      
     size_t x=0;
     for(size_t i=0 ; i<N; i++) 
     {
        x=indices_chosen[i];
        eta[x]=1;
       }
 
 }


template<size_t L> 
void LightLadder<L>::R4WhoActiveIC(const NbrsTable& nbrs)
{
 //cout<<"Calling .R4WhoActiveIC() method: \n" <<endl;   
 NACT=0; //<---SOMETHING BOTHERS ME ABOUT THIS... NACT MEMBER TO CLASS, PRIVATE
 size_t right_nbr=0;
 size_t left_nbr=0;
 size_t across_nbr=0;
 int right_empty=0;
 
 
 for(size_t site=0; site<2*L; site++)
 {
     if(eta[site] == 1)
     {
         right_nbr=nbrs.elems[site][1];
         across_nbr=nbrs.elems[site][2];
         right_empty=2-eta[right_nbr]-eta[across_nbr];
         
         if(right_empty>0)
         {
             left_nbr=nbrs.elems[site][0];
             if(eta[left_nbr]==1 || eta[across_nbr]==1)
             {
                 active_sts[NACT]=site;
                 pos[site]=NACT;
                 NACT++;
                 }
             }
             
         }
          
     }
    
}
/*======================================================================================
 * 
 * Feb 28, 2023, Based on rule4critical.h programmed in 2022. The local variables I 
 * keep reassigning in the OneMove() function make me want to just have a class
 * 
 * RULE 4 ONE MOVE: WE UPDATE SITE SELECTED, BUT THE SITE SELECTED IS NOT SELECTED
   BY THE FOLLOWING METHOD, BUT RATHER JUST ASSIGNED TO LightLadder object, say, we   
   DECLARE LightLadder<L> Lattice and then 
   
    * r_number=RandomSample1(nact, mt) and
    * Lattice.siteSelected= active_sts[r_number]
 *=====================================================================================*/

template<size_t L> 
void LightLadder<L>::R4OneMove(const NbrsTable& nbrs)
{
        
     //cout<<"Calling .R4OneMove() method: \n" <<endl;
    //The followin useful so don't access the O(L) array over and over
     number_locals_changed=0;
     
     old_across_nbr=nbrs.elems[siteSelected][2];
     old_left_nbr=nbrs.elems[siteSelected][0];
     old_right_nbr=nbrs.elems[siteSelected][1];

     
     eta_left=eta[old_left_nbr];
     eta_right=eta[old_right_nbr];
     eta_across=eta[old_across_nbr];
     
      
     
     /* Remember, local_sites_changed[6] is a relative to this class, "global private member" 
      * declared above, maybe has old values, flush them out before*/
 
     for(int i=0; i<6;i++){local_sites_changed[i]=0;}
      

     /*   WE LABEL AND KEEP TRACK OF TYPES OF MOVES in private class member array  
      *   status_updates[] FROM 0-2 AS FOLLOWS
      * 
      *   0. stayed the same
      *   1. went active-->inactive
      *   2. went inactive-->active
      * 
      *   Note since we are dealing with siteSelected, this site just flipped from being
      *   occupied and active, to being empty and inactive, so we update this:
      * */
      
      status_updates[siteSelected]=1; //went active-->inactive
      local_sites_changed[0]=siteSelected;
      number_locals_changed++;
      
      /*-------------------------------------------------------------------------
       * NOW IT SEEMS LIKE I DIVIDED STATUS UPDATE BOOKKEEPING ACCORDING TO
       * WHETHER IT WAS A MOVE TO THE RIGHT, OR AN ACROSS MOVE. SO THE FOLLOWING
       * IF LOOP BEGINS WITH CASE WHERE IT WAS A JUMP TO THE RIGHT 
       *-------------------------------------------------------------------------*/
       
       if(eta_right==0)/*<--note since siteSelected was active, if nbr_right empty,
                             then it will jump to this site now "it prefers move right" */
                             
       {
           
/*-------------------------------------------------------------------------------
            * BEGINNING OF RIGHT MOVE CONDITIONAL IF PART. WHOLE BIG PART.
* ----------------------------------------------------------------------------*/      
     
/* Update info: active_sts, pos, NACT ***/
           
           /*  ---> [FOR MOVE TO THE RIGHT, WE KNOW THAT]: Suppose site 2 is the site chosen, 
            *        site 0 is SITE_LEFT, site 3 is SITE_ACROSS, etc
            *       (for MOVE ACROSS see below ELSE conditional of loop)
            * 
            *        0 2 4 6    <-----(D1) This numbering just for guidance, 1-D type numbering
            *          3 5
            * 
            *          2 4      <-----(D2) This is what numbering looks like if on row below.
            *        1 3 5 7
            * 
Let's see how many right_empty components we can modify. 6 sites involved. * */

           
           eta[old_right_nbr]=1; eta[siteSelected]=0;
           size_t across_col_left=nbrs.elems[old_left_nbr][2];
           int eta_row_across_col_left=eta[across_col_left]; //<--this site plays a "key" role in this case
           // Note because this IF loop only entered if moved to right, then either of old_left or old_across
           // are occupied. In both cases you access site at row_across_col_left (=1 in (D1))
           size_t across_col_right = nbrs.elems[old_across_nbr][1]; 
           int eta_row_across_col_right=eta[across_col_right]; //similar to above
           
           /*=================== [LEFT OLD] ===================*/
           
           if(eta_left==1)
           {
               if(eta_row_across_col_left==1)//<--this would be site# 1 in (D1) above diagram
                   {
                       local_sites_changed[number_locals_changed]=old_left_nbr;
                       number_locals_changed++;
                       status_updates[old_left_nbr]=2;
                       }
               
               }
               
           /*=================== [ACROSS OLD]  ===================*/
           
           if(eta_across==1)
           {
               
               if((eta_row_across_col_left==1) && (eta_row_across_col_right==1))
               {
                   local_sites_changed[number_locals_changed]=old_across_nbr;
                   number_locals_changed++;
                   status_updates[old_across_nbr]=2;
                   }
                else if ((eta_row_across_col_left==0) && (eta_row_across_col_right==0))
                {
                   local_sites_changed[number_locals_changed]=old_across_nbr;
                   number_locals_changed++;
                   status_updates[old_across_nbr]=1;
                    }
  
               }
           
           /*=================== [RIGHT OLD=NEW SITE!]  ===================*/
           
           //By definition eta[RIGHT_OLD]=1
           size_t rightPlus2=nbrs.elems[old_right_nbr][1];
           
           int eta_right_plus2=eta[rightPlus2];
           
           if((eta_row_across_col_right==1)&&(eta_right_plus2==0))
           {
               local_sites_changed[number_locals_changed]=old_right_nbr;
               number_locals_changed++;
               status_updates[old_right_nbr]=2;
               }//<---This is the ONLY TYPE OF UPDATE POSSIBLE!
            
           /*=================== [NEW ACROSS=ROW ACROSS, COL RIGHT]  ===================*/   
            
           if(eta_row_across_col_right==1)
           {
               int eta_right_plus2_across=eta[nbrs.elems[across_col_right][1]];
               if((eta_across==1)&&(eta_right_plus2_across==1))
                   {
                       local_sites_changed[number_locals_changed]=across_col_right;
                       number_locals_changed++;
                       status_updates[across_col_right]=1;
                       }
               else if ((eta_across==0)&&(eta_right_plus2_across==0))
               {
                       local_sites_changed[number_locals_changed]=across_col_right;
                       number_locals_changed++;
                       status_updates[across_col_right]=2;
                   }       
               }

           /*=================== [NEW RIGHT=RIGHT+2]  ===================*/ 

           if(eta_right_plus2==1)
               {
                   int eta_right_plus2_across=eta[nbrs.elems[across_col_right][1]];
                   if(eta_right_plus2_across==0){
                       local_sites_changed[number_locals_changed]=rightPlus2;
                       number_locals_changed++;
                       status_updates[rightPlus2]=2;
                       
                       }
                   }  
           
       }
/*-------------------------------------------------------------------------------
            * END OF RIGHT MOVE CONDITIONAL IF PART.  
* ----------------------------------------------------------------------------*/    

/*******************************************************************************
 *  Here we have 0 2 4 <-----where site 2 is the one that was chosen
 *               1 3 5
 * 
 * Since the new site is now old_across, and the only way to jump across is if
 * eta[old_left]=1=eta[old_right], then we have some info for this update
 ******************************************************************************** */    
       else
       {
           eta[old_across_nbr]=1; eta[siteSelected]=0;
           
           //Useful sites:
           size_t across_col_left=nbrs.elems[old_left_nbr][2];
           int eta_row_across_col_left=eta[across_col_left]; 
           size_t across_col_right = nbrs.elems[old_across_nbr][1]; 
           int eta_row_across_col_right=eta[across_col_right]; //similar to above
           
           /*=================== [LEFT OLD] ===================*/
           
           /*  As said above, eta[left_old]=1 in this across-type move */
           
           if(eta_row_across_col_left==1)
           {
               local_sites_changed[number_locals_changed]=old_left_nbr;
               number_locals_changed++;
               status_updates[old_left_nbr]=2;  
               }
               
           /*=================== [ACROSS OLD=NEW SITE!] ===================*/ 

           if(eta_row_across_col_left==1)
           {
               local_sites_changed[number_locals_changed]=old_across_nbr;
               number_locals_changed++;
               status_updates[old_across_nbr]=2;
               }
           
           /*=================== [RIGHT OLD] ===================*/
           
           if(eta_row_across_col_right==0)
           {
               local_sites_changed[number_locals_changed]=old_right_nbr;
               number_locals_changed++;
               status_updates[old_right_nbr]=1;
               
               }
               
           /*=================== [NEW LEFT=ACROSS LEFT] ===================*/

           if(eta_row_across_col_left==1)//<--automatically went inactive because eta[old_left]=1
            {
                local_sites_changed[number_locals_changed]=across_col_left;
                number_locals_changed++;
                status_updates[across_col_left]=1;   
                }
            
           /*=================== [NEW RIGHT=ACROSS RIGHT] ===================*/
           
            /*!!!! 03/03/2023 UPDATE: NOTHING CHANGES FOR THIS SITE, BECAUSE ETA[RIGHT]=1=ETA[LEFT]
              SINCE THIS WAS AN ACROSS MOVE! Thus new_right, if there is a particle there, it either 
               * remains active, or remains inactive, since a particle at it's left does not affect anything
               * ....
               * 
               * MY OLD PROGRAM HAD AN IF LOOP THAT NOW I NOTICED IS UNNECESSARY*/
    
           
       }
       
       
 
    /*************************  STEP 6: IMPLEMENT DECISIONS: TRICKY PART  **********************
     ******************** UPDATE POS, ACTIVE_STS BUT ONLY UNDER THIS IF-LOOP RIGHT-MOVE  ***************/
            
            /*---->NOTE THIS IS A 2-D version of what previous (to July 2022) Speer-inspired programs do.
             *     we have active_sts[2*L][2], pos[2][L+2], where we only care about active_sites[q][] 
             *     for 0 <= q <= NACT-1 
             * 
             * thus note (active_sts[q][0], active_sts[q][1])=(row of qth active, col of qth active) AND
             * pos[row][col]=q tells you the qth row of active_sts[q][] that you occupy
             * 

 * UPDATES 1 MAIN GIST HERE---> IF A SITE [row_k][col_k] went ACTIVE-->INACTIVE, WE MUST REMOVE IT    
 *                              FROM THE LIST ACTIVE_STS[q][]. We know that a given [row_k][col_k]
                                (in "global" coordinates, not "local" like older versions of program)
 *                              occupies row number q=POS[row_k][col_k] of the array ACTIVE_STS[q][]. 
 * 
 * KEY POINT!!!  ------------->  We REMOVE TYPE 1 UPDATES, AND REPLACE THEM IN THAT SAME ROW q WITH WHAT 
 *                               USED TO BE THE LAST ELEMENT OF ACTIVE_STS[][] THAT ACTUALLY MATTERED,
 *                               NAMELY, ACTIVE_STS[NACT_PREV-1][].
 *                                 
 * 
 * POTENTIAL NOTATION CONFUSION: suppose nact:=nact_previous. Then nact_new=nact-1. 
 *                               Thus active_sts[nact_new][]=previous last "important element" of
 *                               array active_sts[][], since remember it has active_sts[2L][2] dimensions 
 *                               but only the ones from rows in {0,..,q,.,nact-1} matter. 
 *                               
 * 
 * 
 *               In the end, what matters is, if NACT was updated to NACT=NACT-1, 
 *               then old_last_site= (m_row, m_col) 
 *               with m_row=active_sts[nact][0] and m_col=active_sts[nact][1].
 *               With q=pos[row_k][col_k] with (row_k,col_k) a site to remove, we now update  
 *               active_sts to active_sts[q][0]=m_row, active_sts[q][1]=m_col. 
 *               We update pos[m_row][m_col=q to account for this
            
             * */
            
            //SOME DUMMY VARIABLES FOR THIS:
            size_t site_k=0;
            size_t site_m=0;
            size_t q {0}; //<---STAYS SAME DIMENSION as prev to July 2022 programs                
     
 
       for(int k=0; k<number_locals_changed; k++)
       {
           site_k=local_sites_changed[k];
           
           if(status_updates[site_k]==1) //went active-->inactive... thus remove
           {
               NACT--;
               site_m=active_sts[NACT];
               q=pos[site_k];
               active_sts[q]=site_m;
               pos[site_m]=q;
               }
           else //means type 2 update then, since we're only among indices that changed status
           {
               NACT++;
               active_sts[NACT-1]=site_k;
               pos[site_k]=NACT-1;
   
               }    

           }

/*********************************************************
 * END OF UPDATE FOR LOOP, ALSO END OF R4OneMove()
 **********************************************************/
}


/*********************************************************************************
 * *******************************************************************************
 * Now for ONE IC method!! Calls OneMOve() maybe 10^12 times
 * 
 * STARTED 03/03/2023: the below comments are from previous 2022 editions 
 *********************************************************************************
 * ********************************************************************************/


/**************************************************************************************************
 * AUGUST 11, 2022: OneICRule4H2() created to RETURN a harmonic means vector<double> with objectives/characteristics:
 * 
 *  0. BE LIGHT ON MEMORY, as point 2 below details. 

 *  1. We return a vector<double> because it is easy to integrate this vector container within a 
 *     larger array HARMONIC_MEANS[] defined in calling function calling function,
       which is typically a main():
 
 vector<vector<double>> HARMONIC_MEANS(NT, vector<double>(intervals));
   
  * 2.  This function should be lighter than what previous 2021 versions did. In summary we 
  *     eliminate returning or keeping track of a huge rho_a vector which could easily
        have O(10^9) (potentially more) entries. Therefore, inspired on the 
        previously used header "harmonic_mean_ladder.h" function Harmonic_Means2k(), we compute
        and keep track of running sum = \sum_{s=lower}^{upper} [\rho_a(eta_s)]^{-1}
        over intervals type [lower,upper] defined by powers 2^k. But we do not keep the
        individual rho_a(eta_s) terms which become too much. 
           
         *  COMMENT: PREVIOUS 2021 HARMONIC MEANS COMPUTING FUNCTIONS: waited for a whole rho_a
            vector of ULS entries to return, then the main() calling function processed
         *  these entries by calling another function Harmonic_Means2k() to get
         *  a harmonic means vector associated to this eta_0 trial run. 
         * 
 ------------------------------------------------------------------------------------------------------------------------------
  
 Note you have to provide parameter power2 (=k+1old versions) which gives an upper bound on powers
 of 2 and ULS=2^{power2}-1. Kind of confusing but it is for historical reasons I use k=power2-1
  : Previous functions said
  * 
   ULS=2^{k+1}-1, for k an integer. VERY IMPORTANT THAT ULS=2^{k+1}-1 for this to work.
 *  For instance, if k=3, you would have the 2^4-1=15 values of rho_a with indices
 * divided into k+1 containers:
 *                
 *                      || 0 || 1 2 || 3 4 5 6 || 7 8 9 10 11 12 13 14 ||
 *   where the jth (j=0,1,...,k) container has 2^j elements. 

       Note that ULS-1=14=2^4-2, not 2^4-1. This will give k+1 harmonic means
        (the first one,1/( 1/rho_a(0))=rho_a(0) is trivial
 
 
 ************************************************************************************************ */


template<size_t L>
vector<double> LightLadder<L>::OneICR4H2(const int power2, const NbrsTable& nbrs,  mt19937& mt)
{
    //cout<<"Calling .OneICR4H2() method: \n" <<endl;
    int size_means=power2;
    //long double ULS=pow(2, power2)-1;
    vector<double> harmonic_means(size_means);
    long double lower=1;
    long double upper=2;
    long double to_add=2; 
    double running_sum_inverses=0;
    double tempry=0;
    
    harmonic_means[0]=NACT/float(2*L); //<---Note we know NACT member    

/*////////////////////////////////////////////// ///////////////////////////////////////////
                     
                                            DYNAMICAL RULE!
      
///////////////////////////////////////////////////////////////////////////////////////////*/
      
    size_t i; //<---to be the index
      
      
    long double count=0;
      
    if(NACT==0){
          cout<<"Already frozen!"<<'\n';
          }
    else{      
      
      int j=1;
      while(NACT>0&& (j<size_means)) //<---HAD ++j PREVIOUS VERSIONS. 
      {
        
        running_sum_inverses=0;
        //cout<<"note that lower=" <<lower<<" while upper= "<<upper<<endl;
        //cout<<"We have (lower, upper)=("<<lower<<","<< upper<<")"<<'\n';
        while(NACT>0 &&(count<upper))  
        {
                i=RandomSample1(NACT, mt);
                siteSelected=active_sts[i]; //always an active particle selected
                R4OneMove(nbrs); //<--this modifies nact
                count++;
                //cout<< "At step number "<< count << ", we now have NACT=" <<GetNact()<< "\n";
                //PrintEta();
                if(NACT!=0)
                    { 
                      //cout<<"(a) We have nact= "<<nact<<'\n'; 
                      tempry= 2*L/double(NACT);
                      //cout<<"(b) We add= "<<tempry<<'\n';
                      running_sum_inverses+=tempry;
                      //cout<<"(c) atm count="<<count<<'\n';
                    }
                  //else{cout<<"FROZE!"<<'\n';}  
             
            }
            
           harmonic_means[j]=to_add/running_sum_inverses; 
           lower=upper+1;
           to_add=2*to_add;
           upper+=to_add;
           j++;
             
        }
      

      steps_taken=count; 
    }
      //if(NACT==0){cout<< "We FROZE!"<<endl;}
      //cout<< "count="<<count<<endl;
      return harmonic_means;



    
}













 
 

#endif


