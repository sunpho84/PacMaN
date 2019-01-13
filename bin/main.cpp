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
  const vector<array<int,2>> nFreeLegsWhenAssigning;
  
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
  vector<array<int,2>> getNFreeLegsWhenAssigning() const
  {
    /// Returned list
    vector<array<int,2>> out(ass.size());
    
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
  
  /// Gets the list of non-null associations
  vector<NnAss> getNonNullAssociations() const
  {
    /// Result
    vector<NnAss> res;
    
    /// Index of association
    int iAss=
      0;
    
    for(int i=0;i<nPoints;i++)
      for(int j=i+1;j<nPoints;j++)
	{
	  if(ass[iAss]!=0)
	    res.push_back({{i,j},iAss,ass[iAss],nFreeLegsWhenAssigning[iAss]});
	  
	  iAss++;
	}
    
    return
      res;
  }
  
  /// Loops on all Wick contractions, executing the function on it
  template <typename F>
  void forAllWicks(F f) const
  {
    /// Non-null associations
    vector<NnAss> nnAss=
      getNonNullAssociations();
    
    /// Precomputed list of all possible associations
    vector<vector<vector<int>>> possTable(2*nnAss.size());
    
    /// Tensor product of all assignemnt heads and tail case
    vector<int64_t> curr(2*nnAss.size());
    
    /// Last assignemnt before overflow
    vector<int64_t> last(2*nnAss.size());
    
    for(int iNnAss=0;iNnAss<(int)nnAss.size();iNnAss++)
      {
	/// Nonnull ass
	const NnAss& a=
	  nnAss[iNnAss];
	
	/// Index of the association
	const int& iAss=
	  a.iAss;
	
	/// Number of legs to assign
	const int& nLegsToAss=
	  a.nLines;
	
	/// Number of free legs when assignigning the head
	const int& nFreeLegsFrom=
	  nFreeLegsWhenAssigning[iAss][FROM];
	
	/// Last possible assignment on the head
	int64_t lastPossFrom=
	  lastCombination(nLegsToAss,nFreeLegsFrom);
	
	/// Store the starting side combination
	for(int64_t possFrom=firstCombination(nLegsToAss,nFreeLegsFrom);possFrom<lastPossFrom;possFrom=nextCombination(possFrom))
	  possTable[2*iNnAss+FROM].push_back(decryptCombination(nLegsToAss,nFreeLegsFrom,possFrom));
	
	/// Number of free legs when assigning the tail
	const int& nFreeLegsTo=
	  nFreeLegsWhenAssigning[iAss][TO];
	
	/// Last possible assignment on the tail
	int64_t lastPossTo=
	  lastDisposition(nLegsToAss,nFreeLegsTo);
	
	/// Store the ending side disposition
	for(int64_t possTo=firstDisposition(nLegsToAss,nFreeLegsTo);possTo<lastPossTo;possTo=nextDisposition(possTo))
	  possTable[2*iNnAss+TO].push_back(decryptDisposition(nLegsToAss,nFreeLegsTo,possTo));
      }
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
  
  wicksFinder.forAllWicks([](){});
  
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
