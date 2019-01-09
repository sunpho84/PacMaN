#include "Assignment.hpp"

void drawAllAssignments(ostream& os,const vector<Assignment>& all,const vector<int>& nPoint)
{
  /// Gets back N
  const int N=
    nPoint.size();
  
  // Header
  os<<"\\documentclass{standalone}"<<endl;
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
	[iG](int i)
	{
	  return "N_"+to_string(iG)+"_"+to_string(i);
	};
      
      os<<"\\digraph{G"<<iG<<"}{ "<<endl;
      
      // Draws the labels
      for(int iNode=0;iNode<N;iNode++)
	os<<node(iNode)<<" ["
	  <<" label=\""<<nPoint[iNode]<<"\""
	  <<" pos=\""<<cos(2*M_PI*iNode/N)<<","<<sin(2*M_PI*iNode/N)<<"!\""
	  <<" shape=\"circle\""
	  <<"];"<<endl;
      
      // Loop on all assignments
      int i=0;
      for(int row=0;row<N;row++)
	for(int col=row+1;col<N;col++)
	  for(int h=a[i++];h>0;h--)
	    os<<node(row)<<" -> "<<node(col)<<" [arrowhead=none]"<<endl;
      
      // Trailer
      os<<"}"<<endl;
    }
  
  os<<"\\end{document}"<<endl;
}
