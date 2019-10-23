#include "Combinatorial.hpp"

/// List all partitioning of the number m
vector<vector<int>> listAllPartitioningOf(int m)
{
  /// List to be returned
  vector<vector<int>> out;
  
  /// Position
  int k=0;
  
  /// Current vector
  vector<int> p(m,0);
  p[0]=m;
  
  do
    {
      /// Whether to include or not current partition
      bool add=
	true;
      
      // Condition: no 1 must be present
      for(auto i=p.begin();i!=p.begin()+k+1;i++)
      	if(*i==1)
      	  add=false;
      
      // Add or not
      if(add)
	out.emplace_back(p.begin(),p.begin()+k+1);
      
      /// Value to remove
      int remVal=
	0;
      
      // Collect
      while(k>=0 and p[k]==1)
	remVal+=p[k--];
      
      if(k>=0)
	{
	  // Distribute 1
	  p[k]--;
	  remVal++;
	  
	  // Distribute the others
	  while(remVal>p[k])
	    {
	      p[k+1]=p[k];
	      remVal-=p[k++];
	    }
	  
	  p[++k]=remVal;
	}
    }
  while(k>=0);
  
  return out;
}
