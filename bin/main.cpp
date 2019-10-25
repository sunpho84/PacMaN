#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include "Assignment.hpp"
#include "Combinatorial.hpp"
#include "Tools.hpp"
#include "Wick.hpp"

#include <cstring>
#include <fstream>
#include <map>
#include <sstream>

#include <mpi.h>

ofstream realCout("/dev/stdout");
ofstream fakeCout("/dev/null");

/// Number of ranks
int nRanks;

/// Rank id
int rankId;

/// Type used to represent the leg
using S=
  int;

/// Partition of all points, representing a multitrace
vector<Partition<S>> getTraceFromInput(int narg,char **arg)
{
  /// Result
  vector<Partition<S>> pointsTraces(1);
  
  for(int iArg=1;iArg<narg;iArg++)
    {
      /// Temporary scan result
      int t;
      
      /// Result of scanning
      int rc=
	sscanf(arg[iArg],"%d",&t);
      
      bool failed=
	false;
      
      if(rc==1)
	{
	  pointsTraces.back().push_back(t);
	  ostringstream os;
	  os<<t;
	  if(strcasecmp(os.str().c_str(),arg[iArg])!=0)
	    failed=
	      true;
	}
      else
	if(strcasecmp(arg[iArg],",")==0)
	  pointsTraces.emplace_back();
	else
	  failed=
	    true;
      
      if(failed)
	{
	  if(rankId==0)
	    cerr<<"Error! Invalid "<<iArg<<"-th argument "<<arg[iArg]<<", use e.g. "<<arg[0]<<" 2 , 3 , 3"<<endl;
	  MPI_Abort(MPI_COMM_WORLD,0);
	}
    }
  
  COUT<<"Parsed Trace: "<<pointsTraces<<endl;
  
  return
    pointsTraces;
}

/// Count the number of closed loops of the permutation g
inline S countNClosedLoops(vector<S> g)
{
  /// Number of closed loops found
  S nClosedLoops=
    0;
  
  /// Current position
  S i=
    0;
  
  while(i<(int)g.size())
    if(g[i]<0)
      i++;
    else
      {
	S next=g[i];
	// COUT<<i<<endl;
	g[i]=-1;
	i=next;
	
	// COUT<<"Closed loop"<<endl;
	nClosedLoops+=(g[next]<0);
      }
  
  return
    nClosedLoops;
}

/// Transform the points partition into a list of Wick contractions
template <typename S>
Wick<S> makeWickOfPartitions(const vector<Partition<S>>& pointsPart)
{
  /// Result
  Wick<S> out;
  
  /// Running leg
  S iLeg=
    0;
  
  for(auto& pointPart : pointsPart)
    for(auto& part : pointPart)
      {
	for(S i=0;i<part;i++)
	  out.push_back({(S)(iLeg+i),(S)(iLeg+(i+1)%part)});
	iLeg+=part;
      }
  
  return out;
}

/// Returns the string representing the trace part
template <typename S>
string traceDot(const vector<Partition<S>>& pointsTraces)
{
  ostringstream traceNodes;
  
  int u=0;
  int iPT=0;
  
  for(auto pointTrace : pointsTraces)
    {
      traceNodes<<"subgraph cluster_"<<iPT++<<"{"<<endl;
      for(auto p : pointTrace)
	{
	  for(S i=0;i<p;i++)
	    traceNodes<<"\t_"<<u+2*i+1<<" -> _"<<u+2*((i+1)%p)<<endl;
	  u+=2*p;
	}
      
      traceNodes<<"}"<<endl;
    }
  
  return
    traceNodes.str();
}

/// Compute the color factor of this diagram and trace
template <typename S>
void getColFact(S& sign,S& nPow,const S& nLines,const Wick<S>& wick,const int64_t& iCD,vector<S>& totPermSingleContr)
{
  // Count the number of disconnected
  S nDiscoTraces=
    0;
  
  for(S iLine=0;iLine<nLines;iLine++)
    {
      /// Line to consider
      const Line<S>& w=
	wick[iLine];
      
      /// Determine whether the iLine bit is 0 (conn) or 1 (disco)
      const bool CD=
	getBit(iCD,iLine);
      
      // Count the number of disconnected traces, which counts (-1/ncol)^ndisco
      nDiscoTraces+=CD;
      
      /// We swap in1 and in0 if CD is 1
      const S ou0=w[FROM]*2;
      const S ou1=w[TO]*2;
      const S in0=w[FROM^CD]*2+1;
      const S in1=w[TO^CD]*2+1;
      
      totPermSingleContr[ou0]=in1;
      totPermSingleContr[ou1]=in0;
    }
  
  /// Determine the number of closed loops, which counts ncol^nloops
  const S nClosedLoops=
    countNClosedLoops(totPermSingleContr);
  
  /// Parity of nDiscoTraces
  const bool parity=
    nDiscoTraces%2;
  
  sign=
    1-parity*2;
  
  nPow=
    nClosedLoops-nDiscoTraces;
}

int main(int narg,char **arg)
{
  MPI_Init(&narg,&arg);
  
  MPI_Comm_size(MPI_COMM_WORLD,&nRanks);
  
  MPI_Comm_rank(MPI_COMM_WORLD,&rankId);
  
  COUT<<"NRanks: "<<nRanks<<endl;
  
  /// Initial time
  const auto absStart=
    takeTime();
  
  /// Partition of all points, representing a multitrace
  vector<Partition<S>> pointsTraces=
    getTraceFromInput(narg,arg);
  COUT<<"Computing Trace: "<<pointsTraces<<endl;
  
  /// Gets all Wick contractions
  Wick<S> traceStructure=
    makeWickOfPartitions(pointsTraces);
  
  /// Defines the N-Point function
  vector<S> nPoints;
  for(auto p : pointsTraces)
    nPoints.emplace_back(accumulate(p.begin(),p.end(),0));
  
  /// Number of all points
  const S nTotPoints=
    accumulate(nPoints.begin(),nPoints.end(),0);
  
  /// Finder of all assignments
  AssignmentsFinder<S> assignmentsFinder(nPoints);
  
  /// All assignments
  vector<Assignment<S>> allAss=
    assignmentsFinder.getAllAssignements();
  
  /// Draw all assignments
  // ofstream assignmentTex("assignments.tex");
  // drawAllAssignments(assignmentTex,allAss,nPoints);
  
  /// Total number of blob-connecting lines
  const S nLines=
    nTotPoints/2;
  COUT<<nLines<<endl;
  
  /// Compute the number of all Wick contractions
  int64_t nWicksTot=
    computeNTotWicks(allAss,nPoints);
  COUT<<"Total number of Wick contractions: "<<nWicksTot<<endl;
  
  /// Number of possible way to connect or disconnect
  const int64_t nCD=
    (1<<nLines);
  COUT<<"Number of traces options per Wick: "<<nCD<<endl;
  
  /// Number of all color traces to be computed
  int64_t nTotColTraces=
    nWicksTot<<nLines;
  COUT<<"Total number of traces: "<<nTotColTraces<<endl;
  
  /// Time between consecutive prints
  const int timeBetweenPrints=
    10;
  
  /// Number of Wick contraction processed so far
  int64_t nWicksDonePastAss=
    0;
  
  // Loop on all propagator assignment
  for(auto& ass : allAss)
    {
      /// Initial time
      const auto assStart=
	takeTime();
      
      /// Time before next output
      int nSecToNextOutput=
	timeBetweenPrints;
      
      COUT<<"/////////////////////////////////////////////////////////////////"<<endl;
      COUT<<ass<<endl;
      
      /// Color factor computed
      map<int64_t,int64_t> colFact;
      
      /// Lister of all Wick contractions
      WicksFinder<S> wicksFinder(nPoints,ass);
      
      /// Number of Wick contraction of this assignment
      const int64_t nWicksOfThisAss=
	wicksFinder.nAllWickContrs();
      
      /// Workload for this rank
      const auto wl=
	getWorkload(nWicksOfThisAss);
      
      for(int iWick=wl.beg;iWick<wl.end;iWick++)
	{
	  /// Lister of all Wick contractions
	  const Wick<S> wick=
	    wicksFinder.get(iWick);
	  
	  /// Total permutation representing trace + Wick contractions
	  vector<S> totPermSingleContr(2*nTotPoints,-1);
	  
	  // Fill the trace part
	  for(auto p : traceStructure)
	    {
	      const S in=p[0]*2+1;
	      const S out=p[1]*2;
	      totPermSingleContr[in]=out;
	    }
	  
	  // Loop over whether we take connected or disconnected trace for each Wick
	  for(int iCD=0;iCD<nCD;iCD++)
	    {
	      /// Power of the diagram
	      int nPow;
	      
	      /// Sign of the diagram
	      int sign;
	      
	      getColFact(nPow,sign,nLines,wick,iCD,totPermSingleContr);
	      
	      colFact[nPow]+=
		sign;
	    }
	  
	  const auto now=
	    takeTime();
	  
	  const double elapsed=
	    durationInSec(now-assStart);
	  
	  if(elapsed>=nSecToNextOutput)
	    {
	      nSecToNextOutput+=
		timeBetweenPrints;
	      
	      const int64_t nWicksDoneInThisAss=
		(iWick-wl.beg+1)*nRanks;
	      
	      const int64_t nWicksDoneIncludingThisAss=
		nWicksDonePastAss+nWicksDoneInThisAss;
	      
	      const int64_t nWicksResidueOfThisAss=
		nWicksOfThisAss-nWicksDoneInThisAss;
	      
	      const int64_t nWicksResidueTot=
		nWicksTot-nWicksDoneIncludingThisAss;
	      
	      const double timePerWick=
		elapsed/nWicksDoneInThisAss;
	      
	      double timeToEnd=
		nWicksResidueTot*timePerWick;
	      
	      
	      int iQ=0;
	      vector<pair<int,char>> Q{{60,'s'},{60,'m'},{24,'h'},{30,'d'},{12,'M'},{1,'y'}};
	      while(timeToEnd>10 and iQ<(int)Q.size()-1)
		timeToEnd/=Q[iQ++].first;
	      
	      COUT<<
		"NWick done: "<<nWicksDoneInThisAss<<"/"<<nWicksOfThisAss<<", "
		"elapsed time: "<<int(elapsed)<<" s , "
		"expected for this ass: "<<nWicksOfThisAss*timePerWick<<" s , "
		"time to end of this ass: "<<nWicksResidueOfThisAss*timePerWick<<" s, "
		"in total: "<<nWicksTot*timePerWick<<" s , "
		"time to end: "<<timeToEnd<<" "<<Q[iQ].second<<endl;
	    }
	}
      MPI_Barrier(MPI_COMM_WORLD);
      
      // printf("%d done %ld Wick contr\n",omp_get_thread_num(),nDonePerThread);
      
      const auto befRed=
	takeTime();
      COUT<<"Time needed before reduction: "<<durationInSec(befRed-assStart)<<" s"<<endl;
      
      /// Reduce the colFact
      colFact=
	allReduceMap(colFact);
      
      COUT<<"Time needed to redcue: "<<durationInSec(takeTime()-befRed)<<" s"<<endl;
      
      if(rankId==0)
	{
	  for(auto cf : colFact)
	    printf("%+ld*n^(%ld) ",cf.second,cf.first);
	  printf("\n");
	}
    }
  
  // for(int i=0;i<10;i++)
  //   {
  //     COUT<<"/////////////////////////////////////////////////////////////////"<<endl;
  //     COUT<<wicksFinder.decryptWick(i)<<endl;
  //   }
  
  // vector<int> perm{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
  //  for(int iperm=0;iperm<814000;iperm++)
  //    {
  //      int beg=rand()%20;
  //      int end=rand()%(20-beg)+beg;
       
  //      //COUT<<beg<<" "<<end<<endl;
  //      next_permutation(perm.begin()+beg,perm.begin()+end);
  //    }
   
  // ofstream out_perm("/tmp/perm");
  // out_perm<<"digraph G {"<<endl;
  // for(int i=0;i<(int)perm.size();i++)
  //   {
  //     out_perm<<i<<" -> "<<perm[i]<<endl;
  //     COUT<<" "<<perm[i]<<endl;
  //   }
  
  // out_perm<<"}"<<endl;
  
  MPI_Barrier(MPI_COMM_WORLD);
  COUT<<"Total time needed: "<<durationInSec(takeTime()-absStart)<<" s"<<endl;
  
  MPI_Finalize();
  
  return 0;
}
