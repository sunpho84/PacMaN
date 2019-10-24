#include <Tools.hpp>
#include <Wick.hpp>

int64_t computeNTotWicks(const vector<Assignment>& allAss,const vector<int>& nPoints,const bool verbose)
{
  /// Result returned
  int64_t nTotWicks=
    0;
  
  if(verbose)
    cout<<"List of all assignments: "<<endl;
  
  for(auto& ass : allAss)
    {
      const int64_t nWicks=
	WicksFinder(nPoints,ass).nAllWickContrs(false);
      
      if(verbose)
	COUT<<" "<<ass<<" nWick: "<<nWicks<<endl;
      
      nTotWicks+=
	nWicks;
    }
  
  return
    nTotWicks;
}
