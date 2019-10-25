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
template <typename S>
using Assignment=
  vector<S>;

/// Labels the two components of a pair
enum{FROM,TO};

/// Non null assignment
template <typename S>
class NnAss
{
public:
  
  /// Index of the points where the assignment start and ends
  const array<S,2> iPoint;
  
  /// Index of the assignment
  const int iAss;
  
  /// Number of lines
  const S nLines;
  
  /// Number of legs free to assign in start and end
  const array<S,2> nFreeLegsWhenAssigning;
  
  /// Number of possible assignments of the source and sink of the line
  const array<int64_t,2> nPoss;
  
  NnAss(const array<S,2>& iPoint,const int& iAss,const S& nLines,const array<S,2>& nFreeLegsWhenAssigning) :
    iPoint(iPoint),iAss(iAss),nLines(nLines),nFreeLegsWhenAssigning(nFreeLegsWhenAssigning),
    nPoss({nCombinations(nLines,nFreeLegsWhenAssigning[FROM]),nDispositions(nLines,nFreeLegsWhenAssigning[TO])})
  {
  };
  
};

/// Return the index of row, col in the upper triangular part of a marix of size n
template <typename S>
inline int triId(const S& row,const S& col,const S& n)
{
  return n*row-(row+1)*row/2+col-(row+1);
}

/// Assignment of the lines
///
/// Creates all assignments between points
template <typename S>
class AssignmentsFinder
{
  /// N-point function
  vector<S> N;
  
  /// Size of N
  const S nN=
    N.size();
  
  /// Number of assignment to be done
  ///
  /// Number of entries of the upper triangular part of the N X N
  /// matrix
  const int nD=
    (nN-1)*nN/2;
  
public:
  
  AssignmentsFinder(const vector<S>& N) : N(N) {}
  
  /// Iteratively loop on all (row,col) lines, producing all possible assignments
  void getAllAssignments(vector<Assignment<S>>& allAss,Assignment<S>& ass,S row=0,S col=1)
  {
    /// Gets the index of the assignment line
    const int i=
      triId(row,col,nN);
    
    // If we are not looping on the last element
    if(i<nD)
      {
	/// N to be assigned in the starting point
	const S toBeAssFrom=
	  N[row];
	
	/// N to be assigned in the ending point
	const S toBeAssTo=
	  N[col];
	
	/// Number of points in row which would remain unassigned even
	/// if we assigned to it all other points (excluded col)
	const S unassignedInRow=
	  accumulate(N.begin()+col+1,N.end(),toBeAssFrom,minus<int>());
	
	/// Number of points in col which would remain unassigned even
	/// if we assigned to it all other points (excluded row)
	const S unassignedInCol=
	  accumulate(N.begin()+row+1,N.end(),2*toBeAssTo,minus<S>());
	
	/// Minimal value to assign to this (row,col) line
	///
	/// This is equal to the number of lines to be assigned from
	/// the row point, subctracted of the maximum number of lines
	/// that can be assigned to all points following col. If this
	/// number is negative, we have no lower limit on the number
	/// of lines to be assigned.
	const S minAss=
	  max((S)0,
	      max(unassignedInCol,
		  unassignedInRow));
	
	/// Maximal value to assign to this (row,col)
	///
	/// This is the minimum between the number of free points at
	/// start, and end of the line
	const S maxAss=
	  min(toBeAssFrom,toBeAssTo);
	
	// indent(cout,i)<<"(";
	// for(int r=0;r<nN;r++)
	//   {
	//     if(r==row) cout<<RED;
	//     if(r==col) cout<<GREEN;
	//     cout<<N[r];
	//     if(r==row or r==col) cout<<DEFAULT;
	//     if(r!=nN-1) cout<<",";
	//   }
	// cout<<"), Row: "<<row<<", toBeAss(from,to)): ("<<toBeAssFrom<<","<<toBeAssTo<<"), Col: "<<col<<", iD: "<<i<<", minAss: "<<minAss<<", maxAss: "<<maxAss<<endl;
	
	// Increment the column
	S nextCol=
	  (col==nN-1)?(row+2):(col+1);
	
	/// Increment the row
	S nextRow=
	  (col==nN-1)?(row+1):row;
	
	// Loop on all assignment to this element
	for(int a=minAss;a<=maxAss;a++)
	  {
	    // indent(cout,i)<<" Assigning: "<<a<<endl;
	    
	    // assign
	    ass[i]=a;
	    
	    // Temporarily decrease to be assigned
	    N[row]-=a;
	    N[col]-=a;
	    
	    getAllAssignments(allAss,ass,nextRow,nextCol);
	    
	    // Restore to be assigned
	    N[row]+=a;
	    N[col]+=a;
	    
	    // indent(cout,i)<<" Returning: "<<N<<endl;
	  }
      }
    else
      {
	const bool acceptable=
	  N.back()==0;
	
	if(acceptable)
	  {
	    // indent(cout,i)<<" Found good assignment: "<<ass<<endl;
	    // for(int r=0;r<nN;r++)
	    //   {
	    // 	indent(cout,i)<<"  ";
	    // 	for(int c=0;c<r;c++)
	    // 	  cout<<ass[triId(c,r,nN)]<<" ";
	    // 	cout<<"0 ";
	    // 	for(int c=r+1;c<nN;c++)
	    // 	  cout<<ass[triId(r,c,nN)]<<" ";
	    // 	cout<<endl;
	    //   }
	    
	    allAss.push_back(ass);
	  }
	// else
	//   indent(cout,i)<<" Not a good assignment, remained "<<N.back()<<" to be assigned"<<endl;
      }
  }
  
  /// Gets all assignments
  vector<vector<S>> getAllAssignements()
  {
    /// List of all assignments
    vector<Assignment<S>> allAss;
    
    /// List of all assignments
    Assignment<S> ass(nD);
    
    getAllAssignments(allAss,ass);
    
    return
      allAss;
  }
};

template <typename S>
void drawAllAssignments(ostream& os,const vector<Assignment<S>>& all,const vector<S>& nPoint)
{
  /// Gets back N
  const S N=
    nPoint.size();
  
  // Header
  os<<"\\documentclass{standalone}"<<endl;
  os<<"\\usepackage{morewrites}"<<endl;
  os<<"\\usepackage[pdf]{graphviz}"<<endl;
  os<<endl;
  os<<"\\begin{document}"<<endl;
  os<<endl;
  
  for(int iG=0;iG<(int)all.size();iG++)
    {
      /// Assignment
      auto &a=
	all[iG];
      
      auto node=
	[iG](S i)
	{
	  return "N_"+to_string(iG)+"_"+to_string(i);
	};
      
      os<<"\\neatograph{G"<<iG<<"}{ "<<endl;
      
      // Draws the labels
      for(S iNode=0;iNode<N;iNode++)
	os<<node(iNode)<<" ["
	  <<" label=\""<<(char)('a'+iNode)<<"\""// nPoint[iNode]<<"\""
	  <<" pos=\""<<cos(M_PI*(2*iNode+1)/N)<<","<<sin(M_PI*(2*iNode+1)/N)<<"!\""
	  <<" shape=\"circle\""
	  <<"];"<<endl;
      
      // Loop on all assignments
      int i=0;
      for(S row=0;row<N;row++)
	for(S col=row+1;col<N;col++)
	  for(S h=a[i++];h>0;h--)
	    os<<node(row)<<" -- "<<node(col)<<" [arrowhead=none]"<<endl;
      
      // Trailer
      os<<"}"<<endl;
    }
  
  os<<"\\end{document}"<<endl;
}

#endif
