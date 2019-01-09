#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include "Assignment.hpp"
#include "Tools.hpp"

#include <fstream>

/// Creates all Wick contraction, given a n-point fucntion and an assignment
class WicksFinder
{
  /// Number of legs in each point
  const vector<int> nLegsPerPoint;
  
  /// Number of points
  const int nPoints;
  
  /// Total number of legs
  const int nLegs;
  
  /// Point in which the leg is externally attached
  const vector<int> pointOfLeg;
  
  /// Assignment to be considered
  const Assignment ass;
  
  /// Number of permutations of the legs for each point
  const vector<int64_t> nLegsPermPerPoint;
  
  /// Number of permutation of all lines of the given propagator assignements
  const vector<int64_t> nPermPerAss;
  
  /// Number of independent permutations of the legs for each point
  const vector<int64_t> nIndLegsPermPerPoint;
  
  /// Total number of permutation of all legs of all points
  const int64_t nLegsPermAllPoints;
  
  /// Total number of permutations of all propagator assignment
  const int64_t nPermAllAss;
  
public:
  
  /// Compute the number of all Wick contractions
  int64_t nAllWickContrs()
  {
    cout<<" nLegsPermAllPoints "<<nLegsPermAllPoints<<" , nPermAllAss: "<<nPermAllAss<<endl;
    
    return
      nLegsPermAllPoints/nPermAllAss;
  }
  
  /// Given the index of a Wick contraction, decompose into the lines
  vector<int> decryptWick(int64_t iWick) const
  {
    /// Result
    vector<int> res(nLegs);
    
    // Reincorporate the permutations of assignments
    // iWick*=
    //   nPermAllAss;
    
    for(int iPoint=0;iPoint<nPoints;iPoint++)
      {
	cout<<nLegsPermPerPoint[iPoint]<<endl;
	
	/// Index of the permutation
	int iPermOfPoint=
	  iWick%nLegsPermPerPoint[iPoint];
	
	iWick/=nLegsPermPerPoint[iPoint];
	
	cout<<" iperm: "<<iPermOfPoint<<" , "<<decryptPermutation(nLegsPerPoint[iPoint],iPermOfPoint)<<endl;
      }
    
    return pointOfLeg;
    
    return res;
  }
  
  WicksFinder(const vector<int>& nLegsPerPoint,const Assignment& ass) :
    nLegsPerPoint(nLegsPerPoint),
    nPoints(nLegsPerPoint.size()),
    nLegs(summatorial(nLegsPerPoint)),
    pointOfLeg(fillVector<int>(nLegs,[&](const int ileg)
			       {
				 /// Index of the output point
				 int iPoint
				   =0;
				 
				 /// Sum of all legs of points [0,ipoint)
				 int sumLegs
				   =0;
				 
				 while(ileg>=sumLegs+nLegsPerPoint[iPoint])
				   sumLegs+=nLegsPerPoint[iPoint++];
				 
				 return
				   iPoint;
			       })),
    ass(ass),
    nLegsPermPerPoint(transformVector(nLegsPerPoint,factorial<int64_t>)),
    nPermPerAss(transformVector(ass,factorial<int64_t>)),
    nIndLegsPermPerPoint(fillVector<int64_t>(nLegsPerPoint.size(),[&](const int iPoint)
					     {
					       /// Result initialized to the number of permutations of the leg
					       int64_t out=
						 nLegsPermPerPoint[iPoint];
					       
					       // Divides by the product of the number of permutations of all
					       // assignment which start from this point and arrive to higher points
					       for(int jPoint=iPoint+1;jPoint<nPoints;jPoint++)
						 out/=
						   nPermPerAss[triId(iPoint,jPoint,nPoints)];
					       
					       return out;
					     })),
    nLegsPermAllPoints(productorial(nLegsPermPerPoint)),
    nPermAllAss(productorial(nPermPerAss))
  {
    cout<<" ANNA "<<nLegsPerPoint<<endl;
    cout<<" ANNA "<<ass<<endl;
    cout<<" ANNA "<<nLegsPermPerPoint<<endl;
    cout<<" ANNA "<<nIndLegsPermPerPoint<<endl;
  }
  
};

int main(int narg,char **arg)
{
  /// Defines the N-Point function
  const vector<int> nPoint=
    {3,3,2,4};
  
  /// Finder of all assignments
  AssignmentsFinder assignmentsFinder(nPoint);
  
  /// All assigments
  vector<Assignment> allAss=
    assignmentsFinder.getAllAssignements();
  
  //drawAllAssignments("",allAss,nPoint);
  
  WicksFinder wicksFinder(nPoint,allAss.front());
  
  cout<<"NWick of: "<<allAss.front()<<": "<<wicksFinder.nAllWickContrs()<<endl;

  for(int i=0;i<10;i++)
    {
      cout<<"/////////////////////////////////////////////////////////////////"<<endl;
      cout<<wicksFinder.decryptWick(i)<<endl;
    }
  
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
