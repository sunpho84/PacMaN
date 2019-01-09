#ifndef _ASSIGNMENT_HPP
#define _ASSIGNMENT_HPP

#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

using namespace std;

/// Assignement is a vector of int
using Assignment=
  vector<int>;
  
/// Prints an assignent
template <typename T>
ostream& operator<<(ostream& os,const vector<T>& a)
{
  // Open the bracket
  os<<"(";
  
  // Header
  bool h=
    true;
  
  for(auto& i : a)
    {
      if(h)
	h=false;
      else
	os<<",";
      
      os<<i;
    }
  
  os<<")";
  
  return os;
}

/// Draws an Assignment
void drawAllAssignments(ostream& os,const vector<Assignment>& all,const vector<int>& nPoint);

/// Return the index of row, col in the upper triangular part of a marix of size n
inline int triId(const int row,const int col,const int n)
{
  return n*row-(row+1)*row/2+col-(row+1);
}

/// Assignment of the lines
class AssignmentsFinder
{
  /// N-point function
  vector<int> N;
  
  /// Size of N
  const int nN=
    N.size();
  
  /// Number of assignment to be done: upper triangular part of the matrix
  const int nD=
    (nN-1)*nN/2;
  
public:
  
  AssignmentsFinder(const vector<int>& N) : N(N) {}
  
  /// Iteratively loop on all row and col, producing all possible assignments
  void getAllAssignments(vector<Assignment>& allAss,Assignment& ass,int row=0,int col=1)
  {
    /// Gets the index of assignement
    const int i=
      triId(row,col,nN);
    
    // If we are not looping on the last element
    if(i<nD)
      {
	/// N to be assigned
	const int toBeAss=
	  N[row];
	
	/// Minimal value to assign to this row, col
	const int minAss=
	  max(0,accumulate(N.begin()+col+1,N.end(),toBeAss,minus<int>()));
	
	/// Maximal value to assign to this row, col
	const int maxAss=
	  min(toBeAss,N[col]);
	
	cout<<N<<", Row: "<<row<<", toBeAss: "<<toBeAss<<" Col: "<<col<<", iD: "<<i<<", minAss: "<<minAss<<", maxAss: "<<maxAss<<endl;
	
	// Loop on all assignment to this element
	for(int a=minAss;a<=maxAss;a++)
	  {
	    cout<<" "<<a<<endl;
	    
	    // assign
	    ass[i]=a;
	    
	    // Decrease to be assigned
	    N[row]-=a;
	    N[col]-=a;
	    
	    // Increment the column
	    int nextCol=
	      col+1;
	    
	    /// Keep the row
	    int nextRow=
	      row;
	    
	    /// Check if the row is finished
	    const bool rowFinished=
	      nextCol==nN;
	    
	    // Carriage return
	    if(rowFinished)
	      {
		nextRow++;
		
		nextCol=
		  nextRow+1;
	      }
	      getAllAssignments(allAss,ass,nextRow,nextCol);
	      
	      cout<<" Returning: "<<endl;
	      cout<<N<<endl;
	      N[row]+=a;
	      N[col]+=a;
	      cout<<N<<endl;
	  }
      }
    else
      {
	cout<<" Found good assignment: "<<ass<<endl;
	if(N.back()==0)
	  allAss.push_back(ass);
      }
  }
  
  /// Gets all assignments
  vector<vector<int>> getAllAssignements()
  {
    /// List of all assignments
    vector<Assignment> allAss;
    
    /// List of all assignments
    Assignment ass(nD);
    
    getAllAssignments(allAss,ass);
    
    return allAss;
  }
};


#endif
