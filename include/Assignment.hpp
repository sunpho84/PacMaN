#ifndef _ASSIGNMENT_HPP
#define _ASSIGNMENT_HPP

#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include <Combinatorial.hpp>
#include <Tools.hpp>

using namespace std;

inline ostream& indent(ostream& os,int i)
{
  while(i-->=0)
    os<<"|";
  
  return os;
}

/// Assignment is a vector of int
///
/// It stores the upper triangular part of the N X N symmetric matrix,
/// in which each (i,j) entry represents the number of lines going
/// from point i to point j
using Assignment=
  vector<int>;

/// Labels the two components of a pair
enum{FROM,TO};

/// Non null assignment
class NnAss
{
public:
  
  /// Index of the points where the assignment start and ends
  const array<int,2> iPoint;
  
  /// Index of the assignment
  const int iAss;
  
  /// Number of lines
  const int nLines;
  
  /// Number of legs free to assign in start and end
  const array<int,2> nFreeLegsWhenAssigning;
  
  /// Number of possible assignments of the source and sink of the line
  const array<int64_t,2> nPoss;
  
  NnAss(const array<int,2>& iPoint,const int& iAss,const int& nLines,const array<int,2>& nFreeLegsWhenAssigning) :
    iPoint(iPoint),iAss(iAss),nLines(nLines),nFreeLegsWhenAssigning(nFreeLegsWhenAssigning),
    nPoss({nCombinations(nLines,nFreeLegsWhenAssigning[FROM]),nDispositions(nLines,nFreeLegsWhenAssigning[TO])})
  {
  };
  
};

/// Draws an Assignment
void drawAllAssignments(ostream& os,const vector<Assignment>& all,const vector<int>& nPoint);

/// Return the index of row, col in the upper triangular part of a marix of size n
inline int triId(const int row,const int col,const int n)
{
  return n*row-(row+1)*row/2+col-(row+1);
}

/// Assignment of the lines
///
/// Creates all assignments between points
class AssignmentsFinder
{
  /// N-point function
  vector<int> N;
  
  /// Size of N
  const int nN=
    N.size();
  
  /// Number of assignment to be done
  ///
  /// Number of entries of the upper triangular part of the N X N
  /// matrix
  const int nD=
    (nN-1)*nN/2;
  
public:
  
  AssignmentsFinder(const vector<int>& N) : N(N) {}
  
  /// Iteratively loop on all (row,col) lines, producing all possible assignments
  void getAllAssignments(vector<Assignment>& allAss,Assignment& ass,int row=0,int col=1)
  {
    /// Gets the index of the assignment line
    const int i=
      triId(row,col,nN);
    
    // If we are not looping on the last element
    if(i<nD)
      {
	/// N to be assigned in the starting point
	const int toBeAssFrom=
	  N[row];
	
	/// N to be assigned in the ending point
	const int toBeAssTo=
	  N[col];
	
	/// Number of points in row which would remain unassigned even
	/// if we assigned to it all other points (excluded col)
	const int unassignedInRow=
	  accumulate(N.begin()+col+1,N.end(),toBeAssFrom,minus<int>());
	
	/// Number of points in col which would remain unassigned even
	/// if we assigned to it all other points (excluded row)
	const int unassignedInCol=
	  accumulate(N.begin()+row+1,N.end(),2*toBeAssTo,minus<int>());
	
	/// Minimal value to assign to this (row,col) line
	///
	/// This is equal to the number of lines to be assigned from
	/// the row point, subctracted of the maximum number of lines
	/// that can be assigned to all points following col. If this
	/// number is negative, we have no lower limit on the number
	/// of lines to be assigned.
	const int minAss=
	  max(0,
	      max(unassignedInCol,
		  unassignedInRow));
	
	/// Maximal value to assign to this (row,col)
	///
	/// This is the minimum between the number of free points at
	/// start, and end of the line
	const int maxAss=
	  min(toBeAssFrom,toBeAssTo);
	
	indent(cout,i)<<"(";
	for(int r=0;r<nN;r++)
	  {
	    if(r==row) cout<<RED;
	    if(r==col) cout<<GREEN;
	    cout<<N[r];
	    if(r==row or r==col) cout<<DEFAULT;
	    if(r!=nN-1) cout<<",";
	  }
	cout<<"), Row: "<<row<<", toBeAss(from,to)): ("<<toBeAssFrom<<","<<toBeAssTo<<"), Col: "<<col<<", iD: "<<i<<", minAss: "<<minAss<<", maxAss: "<<maxAss<<endl;
	
	// Increment the column
	int nextCol=(col==nN-1)?(row+2):(col+1);
	
	/// Increment the row
	int nextRow=(col==nN-1)?(row+1):row;
	
	// Loop on all assignment to this element
	for(int a=minAss;a<=maxAss;a++)
	  {
	    indent(cout,i)<<" Assigning: "<<a<<endl;
	    
	    // assign
	    ass[i]=a;
	    
	    // Temporarily decrease to be assigned
	    N[row]-=a;
	    N[col]-=a;
	    
	    getAllAssignments(allAss,ass,nextRow,nextCol);
	    
	    // Restore to be assigned
	    N[row]+=a;
	    N[col]+=a;
	    
	    indent(cout,i)<<" Returning: "<<N<<endl;
	  }
      }
    else
      {
	const bool acceptable=N.back()==0;
	
	if(acceptable)
	  {
	    indent(cout,i)<<" Found good assignment: "<<ass<<endl;
	    for(int r=0;r<nN;r++)
	      {
		indent(cout,i)<<"  ";
		for(int c=0;c<r;c++)
		  cout<<ass[triId(c,r,nN)]<<" ";
		cout<<"0 ";
		for(int c=r+1;c<nN;c++)
		  cout<<ass[triId(r,c,nN)]<<" ";
		cout<<endl;
	      }
	    
	    allAss.push_back(ass);
	  }
	else
	  indent(cout,i)<<" Not a good assignment, remained "<<N.back()<<" to be assigned"<<endl;
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

void drawAllAssignments(ostream& os,const vector<Assignment>& all,const vector<int>& nPoint);

#endif
