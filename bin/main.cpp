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

/// Get iBit bit of i
template <typename I>
bool getBit(const I& i,const int& iBit)
{
  return
    (i>>iBit)&1;
}

int countNClosedLoops(const vector<int>& g)
{
  /// Number of closed loops found
  int nClosedLoops=
    0;
  
  /// Current position
  int i=
    0;
  
  vector<bool> visited(g.size(),false);
  while(i<(int)g.size())
    if(visited[i])
      i++;
    else
      {
	// cout<<i<<endl;
	visited[i]=true;
	i=g[i];
	
	if(visited[i])
	  {
	    // cout<<"Closed loop"<<endl;
	    nClosedLoops++;
	  }
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

int main(int narg,char **arg)
{
  /// Partition of all points, representing a multitrace
  vector<Partition> pointsTraces=
  //{{2},{2},{4},{2,2}};
    {{3},{3},{6},{6}};
  
  Wick traceStructure=
    makeWickOfPartitions(pointsTraces);
  
  ostringstream traceNodes;
  
  {
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
  }
  
  // for(auto y : listAllPartitioningOf(4))
  //   cout<<y<<endl;
  
  // for(auto a : std::vector<std::pair<int,int>>{{4,0},{4,1},{4,2},{4,3},{4,4}})
  //   cout<<newtonBinomial(a.first,a.second)<<endl;
  
  /// Defines the N-Point function
  const vector<int> nPoints=
    //    {2,2,4,4};
     {3,3,6,6};
    // {3,3,3,3};
    //{2,2,2,2,4};
  
  const int nTotPoints=accumulate(nPoints.begin(),nPoints.end(),0);
  
  /// Finder of all assignments
  AssignmentsFinder assignmentsFinder(nPoints);
  
  /// All assigments
  vector<Assignment> allAss=
    assignmentsFinder.getAllAssignements();
  
  ofstream assignmentTex("assignments.tex");
  drawAllAssignments(assignmentTex,allAss,nPoints);
  
  /// Total permutation representing trace + Wick contractions
  vector<int> totPermSingleContr(2*nTotPoints,-1);
  cout<<"Printing the trace"<<endl;
  for(auto p : traceStructure)
    {
      cout<<p[0]*2+1<<" "<<p[1]*2<<endl;
      totPermSingleContr[p[0]*2+1]=p[1]*2;
    }
  
  // Loop on all propagator assignment
  for(auto& ass : allAss)
    {
      cout<<"/////////////////////////////////////////////////////////////////"<<endl;
      cout<<ass<<endl;
      
      /// Lister of all Wick contractions
      WicksFinder wicksFinder(nPoints,ass);
      
      // cout<<"Wick contractions of assignment: "<<ass<< endl;
      
      int64_t iWick=
	0;
      
      map<int,int> colFact;
      
      wicksFinder.forAllWicks([&iWick,&totPermSingleContr,// &traceNodes,
			       &colFact](Wick& wick)
			      {
				//  cout<<endl;
				
				// cout<<"Wick contraction "<<iWick<<endl;
				// for(auto& w : wick)
				//   cout<<" Assigning leg "<<w[FROM]<<" to "<<w[TO]<<endl;
				
				/// Total number of lines
				const int nLines=wick.size();
				
				/// Number of possible way to connect or disconnect
				const int64_t nCD=(1<<nLines);
				
				// Loop over whether we take connected or disconnected trace for each Wick
				for(int iCD=0;iCD<nCD;iCD++)
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
					
					// cout<<"\t _"<<ou0<<" -> _"<<in1<<endl;
					// cout<<"\t _"<<ou1<<" -> _"<<in0<<endl;
				      }
				    
				    /// Determine the number of closed loops, which counts ncol^nloops
				    const int nClosedLoops=
				      countNClosedLoops(totPermSingleContr);
				    
				    //cout<<"Nclosed loops: "<<nClosedLoops<<", ndisco trace: "<<nDiscoTraces<<endl;
				    
				    const bool parity=
				      nDiscoTraces%2;
				    
				    const int nPow=
				      nClosedLoops-nDiscoTraces;
				    //cout<<" "<<(parity?"-":"+")<<"nC^"<<nPow<<endl;
				    
				    colFact[nPow]+=1-parity*2;
				  }
				iWick++;
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
