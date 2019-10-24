#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include "Assignment.hpp"
#include "Combinatorial.hpp"
#include "Tools.hpp"
#include "Wick.hpp"

#include <fstream>
#include <map>
#include <sstream>

#include <omp.h>

/// Get iBit bit of i
template <typename I>
bool getBit(const I& i,const int& iBit)
{
  return
    (i>>iBit)&1;
}

inline int countNClosedLoops(vector<int> g)
{
  /// Number of closed loops found
  int nClosedLoops=
    0;
  
  /// Current position
  int i=
    0;
  
  while(i<(int)g.size())
    if(g[i]<0)
      i++;
    else
      {
	int next=g[i];
	// cout<<i<<endl;
	g[i]=-1;
	i=next;
	
	// cout<<"Closed loop"<<endl;
	nClosedLoops+=(g[next]<0);
      }
  
  return
    nClosedLoops;
}

/// Transform the points partition into a list of Wick contractions
Wick makeWickOfPartitions(const vector<Partition>& pointsPart)
{
  /// Result
  Wick out;
  
  /// Running leg
  int iLeg=
    0;
  
  for(auto& pointPart : pointsPart)
    for(auto& part : pointPart)
      {
	for(int i=0;i<part;i++)
	  out.push_back({iLeg+i,iLeg+(i+1)%part});
	iLeg+=part;
      }
  
  return out;
}

/// Returns the string representing the trace part
string traceDot(const vector<Partition>& pointsTraces)
{
  ostringstream traceNodes;
  
  int u=0;
  int iPT=0;
  
  for(auto pointTrace : pointsTraces)
    {
      traceNodes<<"subgraph cluster_"<<iPT++<<"{"<<endl;
      for(auto p : pointTrace)
	{
	  for(int i=0;i<p;i++)
	    traceNodes<<"\t_"<<u+2*i+1<<" -> _"<<u+2*((i+1)%p)<<endl;
	  u+=2*p;
	}
      
      traceNodes<<"}"<<endl;
    }
  
  return
    traceNodes.str();
}

auto getColFact(const int& nLines,const Wick& wick,const int64_t& iCD,vector<int>& totPermSingleContr)
{
  // cout<<" Writing conn disco "<<iCD<<" = ";
  // for(int iLine=0;iLine<nLines;iLine++)
  //   cout<<((iCD>>iLine)&1);
  // cout<<endl<<endl;
  
  // cout<<traceNodes.str()<<endl;
  
  // Count the number of disconnected
  int nDiscoTraces=
    0;
  
  for(int iLine=0;iLine<nLines;iLine++)
    {
      /// Line to consider
      const Line& w=
	wick[iLine];
      
      /// Determine whether the iLine bit is 0 (conn) or 1 (disco)
      const bool CD=
	getBit(iCD,iLine);
      
      // Count the number of disconnected traces, which counts (-1/ncol)^ndisco
      nDiscoTraces+=CD;
      
      /// We swap in1 and in0 if CD is 1
      const int ou0=w[FROM]*2;
      const int ou1=w[TO]*2;
      const int in0=w[FROM^CD]*2+1;
      const int in1=w[TO^CD]*2+1;
      
      totPermSingleContr[ou0]=in1;
      totPermSingleContr[ou1]=in0;
    }
  
  /// Determine the number of closed loops, which counts ncol^nloops
  const int nClosedLoops=
    countNClosedLoops(totPermSingleContr);
  
  const bool parity=
    nDiscoTraces%2;
  
  const int sign=
    1-parity*2;
  
  const int nPow=
    nClosedLoops-nDiscoTraces;
  
  return
    make_tuple(nPow,sign);
}

int main(int narg,char **arg)
{
  int nThreads;
  
#pragma omp parallel
  nThreads=
    omp_get_num_threads();
  
  /// Partition of all points, representing a multitrace
  vector<Partition> pointsTraces=
    //{{2},{2},{4},{2,2}};
    {{6},{6},{6},{6}};
    // {{3},{3},{6},{6}};
  
  Wick traceStructure=
    makeWickOfPartitions(pointsTraces);
  
  // for(auto y : listAllPartitioningOf(4))
  //   cout<<y<<endl;
  
  // for(auto a : std::vector<std::pair<int,int>>{{4,0},{4,1},{4,2},{4,3},{4,4}})
  //   cout<<newtonBinomial(a.first,a.second)<<endl;
  
  /// Defines the N-Point function
  const vector<int> nPoints=
    //{2,2,4,4};
    {6,6,6,6};
    // {3,3,3,3};
    //{2,2,2,2,4};
  
  /// Number of all points
  const int nTotPoints=
    accumulate(nPoints.begin(),nPoints.end(),0);
  
  /// Finder of all assignments
  AssignmentsFinder assignmentsFinder(nPoints);
  
  /// All assignments
  vector<Assignment> allAss=
    assignmentsFinder.getAllAssignements();
  
  /// Draw all assignments
  ofstream assignmentTex("assignments.tex");
  drawAllAssignments(assignmentTex,allAss,nPoints);
  
  /// Total permutation representing trace + Wick contractions
  vector<int> totPermSingleContr(2*nTotPoints,-1);
  
  // Fill the trace part
  for(auto p : traceStructure)
    {
      const int in=p[0]*2+1;
      const int out=p[1]*2;
      totPermSingleContr[in]=out;
    }
  
  /// Total number of blob-connecting lines
  const int nLines=
    nTotPoints/2;
  cout<<nLines<<endl;
  
  /// Compute the number of all Wick contractions
  int64_t nTotWicks=
    computeNTotWicks(allAss,nPoints);
  cout<<"Total number of Wick contractions: "<<nTotWicks<<endl;
  
  int64_t nTotColTraces=
    nTotWicks<<nLines;
  cout<<"Total number of traces: "<<nTotColTraces<<endl;
  
  // auto wc=
  //   getColFact(nLines,WicksFinder(nPoints,allAss.front()).getFirst(),0,totPermSingleContr);
  // cout<<get<0>(wc)<<" "<<get<1>(wc)<<endl;
  // return 0;
  
  const auto start=
    takeTime();
  
  const int timeBetweenPrints=
    10;
  
  int nSecToNextOutput=
    0;
  
  /// Index od the Wick contraction done so far
  int64_t iWick=
    0;
  
  // Loop on all propagator assignment
  for(auto& ass : allAss)
    {
      cout<<"/////////////////////////////////////////////////////////////////"<<endl;
      cout<<ass<<endl;
      
      /// Lister of all Wick contractions
      WicksFinder wicksFinder(nPoints,ass);
      
      // cout<<"Wick contractions of assignment: "<<ass<< endl;
      
      map<int,int> colFact;
      
      wicksFinder.forAllWicks([&](Wick& wick)
			      {
				//  cout<<endl;
				
				// cout<<"Wick contraction "<<iWick<<endl;
				// for(auto& w : wick)
				//   cout<<" Assigning leg "<<w[FROM]<<" to "<<w[TO]<<endl;
				
				/// Number of possible way to connect or disconnect
				const int64_t nCD=
				  (1<<nLines);
				
				// Loop over whether we take connected or disconnected trace for each Wick
#pragma omp parallel for
				for(int iCD=0;iCD<nCD;iCD++)
				  {
				    int nPow,sign;
				    
				    tie(nPow,sign)=
				      getColFact(nLines,wick,iCD,totPermSingleContr);
				    
				    colFact[nPow]+=
				      sign;
				  }
				
				iWick++;
				if(iWick%10000==0)
				  {
				    const auto now=
				      takeTime();
				    
				    const int nSecFromStart=
				      durationInSec(now-start);
				    
				    if(nSecFromStart>=nSecToNextOutput)
				      {
					nSecToNextOutput=
					  nSecFromStart+timeBetweenPrints;
					
					const double elapsed=
					  durationInSec(now-start);
					
					cout<<
					  "NWick done: "<<iWick<<"/"<<nTotWicks<<", "
					  "elapsed time: "<<elapsed<<" s , "
					  "tot expected: "<<nTotWicks*elapsed/iWick/nThreads<<" s , "
					  "time to end: "<<(nTotWicks-iWick)*elapsed/iWick/nThreads<<" s"<<endl;
				      }
				    
				  }
			      });
      
      for(auto cf : colFact)
	printf("%+d*n^(%d) ",cf.second,cf.first);
      printf("\n");
      
      cout<<"NWick: "<<iWick<<endl;
    }
  
  // for(int i=0;i<10;i++)
  //   {
  //     cout<<"/////////////////////////////////////////////////////////////////"<<endl;
  //     cout<<wicksFinder.decryptWick(i)<<endl;
  //   }
  
  // vector<int> perm{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
  //  for(int iperm=0;iperm<814000;iperm++)
  //    {
  //      int beg=rand()%20;
  //      int end=rand()%(20-beg)+beg;
       
  //      //cout<<beg<<" "<<end<<endl;
  //      next_permutation(perm.begin()+beg,perm.begin()+end);
  //    }
   
  // ofstream out_perm("/tmp/perm");
  // out_perm<<"digraph G {"<<endl;
  // for(int i=0;i<(int)perm.size();i++)
  //   {
  //     out_perm<<i<<" -> "<<perm[i]<<endl;
  //     cout<<" "<<perm[i]<<endl;
  //   }
  
  // out_perm<<"}"<<endl;
  
  return 0;
}
