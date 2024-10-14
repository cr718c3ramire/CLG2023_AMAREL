#ifndef AUXCLG2023_H
#define AUXCLG2023_H
#include<iostream>
#include<cmath>
#include<cstddef>

#include<random>
#include<algorithm>
#include<iterator>

using namespace std;


/********************************************************************************
 * the following function only works for enumeration row-wise order such as in
         1 2 3 4
         5 6 7 8
 *******************************************************************************/
template<typename T, typename T2, size_t L>
size_t index( T x, T2 y )  { return x*L + y; }
//size_t index( size_t x, size_t y )  { return x + L * y; }
//const size_t indexConstants( const int x, const int y )  { return L*x +  y; }
template<typename T>
inline int IsOdd(const T x) { return x & 1; }// uses BITwise AND, faster than division

/********************************************************************************
 * the following function only works for enumeration columnw-wise(MATLAB)
   order such as in
    
         0 2 4 6
         1 3 5 7
 *******************************************************************************/
template<typename T, typename T2>
size_t Index1D( T x, T2 y )  
{
 size_t result=0;
 
 if(x==0){result=2*y;} 
 else{result=2*y+1;}
   
    
 return result; 
 }

/***********************************************************************************
 * The following used to return a vector, but I don't like that anymore. I'll make it
   take 3 arguments, and it has to modify by reference two of the inputs
 * 
 * 
 * *********************************************************************************/ 
template<typename T, typename T1, typename T2>
void Index2D(T x , T1& ex, T2& why)  
{
 
  auto m=x/2;
  if(IsOdd(x)==1)
    {
        ex=1;
        why=floor(m);
        }
  else
    {
        ex=0;
        why=floor(m); 
      }
 }

/*************************************************************************************
 *  SOME FUNCTIONS TO DISPLAY OR PRINT TO SCREEN, RECYCLED FROM OLDER 2021-2022 PROGRAMS
 * * *********************************************************************************/
template<typename OtherVectorClass>

void DisplayVector(OtherVectorClass& vec){
  cout<< "[";  
  for(auto & iterator:vec)
    {
            cout << iterator << ", ";
    }
  cout<< "]"<<'\n';
}

template <typename vec1D, typename file> //WITH OR WITHOUT THIS, could not make following overload
void VectorToFileCPP(vec1D const & vec, file& infile){
  for(auto& iterator: vec)
    {
      infile<< iterator<<"\t"; //double quotes are strings, single quotes are characters. They're differenT! 
    }
  cout<< "\n"<<endl;
}

//This template below assumes you are conforming to STL algorithms/iterators notations
template<typename OtherVectorClass>

void PythonPrintVector(OtherVectorClass& vec){
    int i=0;
    cout<<"np.array([";
  for(auto it = std::begin(vec); it != std::end(vec)-1; ++it)
    {
        
            //cout << iterator <<", ";
            cout << *it <<", ";
            //i++;
    }
  
  cout<< *(std::end(vec)-1) <<"])\n"<<endl;
  
}


template <typename vec1D, typename file>
void PythonFILEPrintVector(vec1D const & vec, file& infile){
    int i=0;
    infile<<"np.array([";
  for(auto it = std::begin(vec); it != std::end(vec)-1; ++it)
    {
        
            //cout << iterator <<", ";
            infile << *it <<", ";
            //i++;
    }
  infile<< *(std::end(vec)-1) <<"])\n"<<endl;
}

size_t RandomSample1(size_t& n, mt19937& mt) //random samples 1 from n, among 0,..,n-1 integers, using seed and mt19937 engine 
// NEED TO input a seed, such as  auto seed=time(0); //sets seed=time on computer
{
     vector<size_t> indices_to_shuffle (n); 
     for (size_t i=0; i<n; ++i){ indices_to_shuffle[i]=i;}
 
      shuffle (indices_to_shuffle.begin(), indices_to_shuffle.end(), mt);
      size_t selected_integer=indices_to_shuffle[0];
      
      return selected_integer;
   
    } 

#endif