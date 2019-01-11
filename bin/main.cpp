#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include "Assignment.hpp"
#include "Combinatorial.hpp"
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
  
  /// Number of permutation of all lines of the given propagator assignments
  const vector<int64_t> nPermPerAss;
  
  /// Number of independent permutations of the legs for each point
  const vector<int64_t> nIndLegsPermPerPoint;
  
  /// Total number of permutation of all legs of all points
  const int64_t nLegsPermAllPoints;
  
  /// Number of free legs in the head and in the tail when assigning the (i,j) assignment
  const vector<pair<int,int>> nFreeLegsWhenAssigning;
  
  /// Total number of permutations of all propagator assignment
  const int64_t nPermAllAss;
  
  /// Compute the number of free legs on the head and tail of all
  /// associations when assigning each of them in turn
  ///
  /// The number of free legs in the head is given by
  /// nLegsPerPoint[row]-\sum_{k=0}^{row-1} ass[k,row]-
  /// \sum_{k=row+1}^{col-1} ass[row,k].
  ///
  /// The number of free legs in the tail is given by
  /// nLegsPerPoint[col]-\sum_{k=0}^{row-1} ass[k,col]
  vector<pair<int,int>> getNFreeLegsWhenAssigning() const
  {
    /// Returned list
    vector<pair<int,int>> out(ass.size());
    
    /// Association considered
    int iAss=
      0;
    
    for(int row=0;row<nPoints;row++)
      for(int col=row+1;col<nPoints;col++)
	{
	  /// Number of free legs in head (row)
	  int nFreeLegsInHead=
	    nLegsPerPoint[row];
	  
	  /// Number of free legs in tail (col)
	  int nFreeLegsInTail=
	    nLegsPerPoint[col];
	  
	  for(int k=0;k<row;k++)
	    {
	      nFreeLegsInHead-=
		ass[triId(k,row,nPoints)];
	      
	      nFreeLegsInTail-=
		ass[triId(k,col,nPoints)];
	    }
	  
	  for(int k=row+1;k<col;k++)
	    nFreeLegsInHead-=
	      ass[triId(row,k,nPoints)];
	  
	  // Set the value
	  out[iAss]=
	    {nFreeLegsInHead,nFreeLegsInTail};
	  
	  iAss++;
	}
    
    return
      out;
  }
  
public:
  
  /// Compute the number of all Wick contractions
  int64_t nAllWickContrs()
  {
    cout<<" nLegsPermAllPoints "<<nLegsPermAllPoints<<" , nPermAllAss: "<<nPermAllAss<<endl;
    
    return
      nLegsPermAllPoints/nPermAllAss;
  }
  
  /// Loops on all Wick contractions, executing the function on it
  void forAllWicks(int64_t iWick) const
  {
    
    
    // /// Result
    // //vector<int64_t> res(nLegs);
    
    // // Reincorporate the permutations of assignments
    // // iWick*=
    // //   nPermAllAss;
    
    // /// Index of the permutation for each point
    // vector<int64_t> iPermLegsOfPoint=
    //   decomposeNumber(iWick,nIndLegsPermPerPoint);
    
    // vector<int64_t> p;
    
    // /// Index running on all legs
    // int iLeg=
    //   0;
    
    // // Loop on all points
    // for(int iPoint=0;iPoint<nPoints;iPoint++)
    //   {
    // 	/// Index running on all legs of this point
    // 	int iLegOfPoint=
    // 	  0;
	
    // 	// Loop on all other points, to run on all associations
    // 	for(int jPoint=0;jPoint<nPoints;jPoint++)
    // 	  if(iPoint!=jPoint)
    // 	    {
    // 	      /// Gets the relevant association
    // 	      const int iAss=
    // 		(iPoint<jPoint)
    // 		?
    // 		triId(iPoint,jPoint,nPoints)
    // 		:
    // 		triId(jPoint,iPoint,nPoints);
	      
    // 	      /// Number of permutations for this side of the
    // 	      /// association. This is equal to the number of
    // 	      /// permutations of the association, or 1 if the lines
    // 	      /// goes out
    // 	      const int64_t nPermLegsOfAss=
    // 		(jPoint>iPoint)
    // 		?
    // 		1
    // 		:
    // 		nPermPerAss[iAss];
	      
    // 	      // Loop on all legs of this association
    // 	      for(int iLegOfAss=0;iLegOfAss<ass[iAss];iLeg++)
    // 		{
    // 		  /// Number of legs to be assigned
    // 		  const int nLegsFree=
    // 		    nLegsPerPoint[iPoint]-iLegOfPoint;
    // 		}
    // 	    }
	
    // //ora tocca decomporre le singole permutazioni nelle permutazioni delle sottoassegnazioni
    // //questo si fa facendo il loop sulle assegnazioni
    
    // // for(int iPoint=0;iPoint<nPoints;iPoint++)
    // //   {
    // // 	cout<<nLegsPermPerPoint[iPoint]<<endl;
	
    // // 	/// Index of the permutation
    // // 	int iPermOfPoint=
    // // 	  iWick%nLegsPermPerPoint[iPoint];
	
    // // 	iWick/=nLegsPermPerPoint[iPoint];
	
    // // 	cout<<" iperm: "<<iPermOfPoint<<" , "<<decryptPermutation(nLegsPerPoint[iPoint],iPermOfPoint)<<endl;
    // //   }
    
    // // return pointOfLeg;
    
    // // return res;
    //   }
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
    nFreeLegsWhenAssigning(getNFreeLegsWhenAssigning()),
    nPermAllAss(productorial(nPermPerAss))
  {
    cout<<" ANNA propStr: "<<nLegsPerPoint<<endl;
    cout<<" ANNA ass: "<<ass<<endl;
    cout<<" ANNA nLegsPermPerPoint: "<<nLegsPermPerPoint<<endl;
    cout<<" ANNA nIndLegsPermPerPoint: "<<nIndLegsPermPerPoint<<endl;
    cout<<" ANNA nFreeLegsWhenAssigning: "<<nFreeLegsWhenAssigning<<endl;
  }
  
};

int main(int narg,char **arg)
{
  forAllCombinations(5, 6, [](int64_t i){cout<<i<<endl;});
  
  return 0;
  
  int nSlots=
    4;
  
  int nObj=
    3;
  
  int nDisp=
    nDispositions(nObj,nSlots);
  
  for(int64_t iDisp=0;iDisp<nDisp;iDisp++)
    cout<<decryptDisposition(nObj, nSlots,iDisp)<<endl;
  return 0;
  
  for(auto a : std::vector<std::pair<int,int>>{{4,0},{4,1},{4,2},{4,3},{4,4}})
    cout<<newtonBinomial(a.first,a.second)<<endl;
  
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
