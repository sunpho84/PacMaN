#ifndef _WICK_HPP
#define _WICK_HPP

#include <memory>

#include "Assignment.hpp"

using namespace std;

/// A line connects two points
template <typename S>
using Line=
  array<S,2>;

/// Wick contraction is a list of lines
template <typename S>
using Wick=vector<Line<S>>;

/// Creates all Wick contraction, given a n-point function and an assignment
template <typename S>
class WicksFinder
{
  /// Number of legs in each point
  const vector<S> nLegsPerPoint;
  
  /// Number of points
  const S nPoints;
  
  /// Total number of legs
  const S nLegs;
  
  /// Total number of lines, equal to half the legs
  const S nLines;
  
  /// Point where the leg is attached
  const vector<S> pointOfLeg;
  
  /// Assignment to be considered
  const Assignment<S> ass;
  
  /// Number of permutations of the legs of each point
  const vector<int64_t> nLegsPermPerPoint;
  
  /// Number of permutations of all lines of the given propagator
  /// assignments, for each point
  const vector<int64_t> nPermPerAss;
  
  /// Total number of permutations of all legs of all points
  const int64_t nLegsPermAllPoints;
  
  /// Number of free legs in the head and in the tail when assigning the (i,j) assignment
  const vector<array<S,2>> nFreeLegsWhenAssigning;
  
  /// Product of the number of all permutations of all lines of the
  /// given propagator assignment
  const int64_t nPermAllAss;
  
  /// Number of legs before the given point
  const vector<S> nLegsBefPoint;
  
  /// Compute the number of free legs on the head and tail of all
  /// assignment when assigning each of them in turn
  ///
  /// The number of free legs in the head is given by
  /// nLegsPerPoint[row]-\sum_{k=0}^{row-1} ass[k,row]-
  /// \sum_{k=row+1}^{col-1} ass[row,k].
  ///
  /// The number of free legs in the tail is given by
  /// nLegsPerPoint[col]-\sum_{k=0}^{row-1} ass[k,col]
  vector<array<S,2>> getNFreeLegsWhenAssigning() const
  {
    /// Returned list
    vector<array<S,2>> out(ass.size());
    
    /// Considered assignment
    int iAss=
      0;
    
    for(S row=0;row<nPoints;row++)
      for(S col=row+1;col<nPoints;col++)
	{
	  /// Number of free legs in head (row)
	  S nFreeLegsInHead=
	    nLegsPerPoint[row];
	  
	  /// Number of free legs in tail (col)
	  S nFreeLegsInTail=
	    nLegsPerPoint[col];
	  
	  for(S k=0;k<row;k++)
	    {
	      nFreeLegsInHead-=
		ass[triId(k,row,nPoints)];
	      
	      nFreeLegsInTail-=
		ass[triId(k,col,nPoints)];
	    }
	  
	  for(S k=row+1;k<col;k++)
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
  
  /// Precomputed list of all possible assignment
  vector<vector<vector<S>>> possTable;
  
  /// Non-null associations
  vector<NnAss<S>> nnAss;
  
  /// Looper on all possibilities
  unique_ptr<Digits<S>> possibilitiesLooper;
  
public:
  
  /// Return first Wick contraction
  Wick<S> getFirst()
  {
    possibilitiesLooper->setTo(0);
    
    return
      convertDigitsToWick(possibilitiesLooper->digits);
  }
  
  /// Compute the number of all Wick contractions
  int64_t nAllWickContrs(const bool verbose=true)
  {
    if(verbose)
      COUT<<" nLegsPermAllPoints "<<nLegsPermAllPoints<<" , nPermAllAss: "<<nPermAllAss<<endl;
    
    return
      nLegsPermAllPoints/nPermAllAss;
  }
  
  /// Gets the list of non-null associations
  vector<NnAss<S>> getNonNullAssociations() const
  {
    /// Result
    vector<NnAss<S>> res;
    
    /// Index of association
    int iAss=
      0;
    
    for(S i=0;i<nPoints;i++)
      for(S j=i+1;j<nPoints;j++)
	{
	  if(ass[iAss]!=0)
	    res.push_back({{i,j},iAss,ass[iAss],nFreeLegsWhenAssigning[iAss]});
	  
	  iAss++;
	}
    
    return
      res;
  }
  
  /// Convert the digits of the Wick contraction id written in terms of digits into an actual Wick contraction
  Wick<S> convertDigitsToWick(const vector<S>& wickDigits)
    const
  {
    /// Store wether the leg is assigned
    vector<int> legIsAss(nLegs,false);
    
    /// Store the assignment of the legs, in form of lines connecting two legs
    Wick<S> lineAss(nLines);
    
    /// Index of the line to assign
    S iLineToAss=
      0;
    
    for(int iNnAss=0;iNnAss<(int)nnAss.size();iNnAss++)
      {
	/// Number of legs for this assignment
	const S nLegsPerAss=
	  nnAss[iNnAss].nLines;
	
	/// Legs connected by the lines of this assignment
	///
	/// We need to temporarily store this, because in the table
	/// of possibility they are counted as a whole block
	vector<array<S,2>> pointsLegAss(nLegsPerAss);
	
	for(S iLine=0;iLine<nLegsPerAss;iLine++)
	  for(int ft=0;ft<2;ft++)
	    {
	      // At first, set the leg to the number of legs preceeding
	      // the point (which is the lable of the first leg of the point)
	      S& l=
		pointsLegAss[iLine][ft]=
		nLegsBefPoint[nnAss[iNnAss].iPoint[ft]];
	      
	      // Then skip needed unassigned legs
	      S count=
		possTable[2*iNnAss+ft][wickDigits[2*iNnAss+ft]][iLine];
	      while(count>0 or legIsAss[l])
		{
		  if(not legIsAss[l])
		    count--;
		  l++;
		}
	    }
	
	// Now we mark all assigned
	for(S iLegPerAss=0;iLegPerAss<nLegsPerAss;iLegPerAss++)
	  {
	    for(int ft=0;ft<2;ft++)
	      {
		lineAss[iLineToAss][ft]=
		  pointsLegAss[iLegPerAss][ft];
		
		legIsAss[pointsLegAss[iLegPerAss][ft]]=
		  true;
	      }
	    
	    iLineToAss++;
	  }
      }
    
    return
      lineAss;
  }
  
  /// Loops on all Wick contractions, executing the function on it
  template <typename F>
  void forAllWicks(F f)
    const
  {
    possibilitiesLooper->forAllNumbers([&,this](const vector<S>& wickDigits)
				       {
					 /// Line assigments
					 Wick<S> lineAss=
					   convertDigitsToWick(wickDigits);
					 
					 f(lineAss);
				       });
  }
  
  /// Get the Wick contraction numberiWick
  Wick<S> get(const int64_t& iWick)
  {
    possibilitiesLooper->setTo(iWick);
    
    return
      convertDigitsToWick(possibilitiesLooper->digits);
  }
  
  /// Reset the WicksFinder
  void reset()
  {
    nnAss=
      getNonNullAssociations();
    
    // for(auto& n : nnAss)
    //   {
    // 	cout<<"Assignment "<<n.iAss<<endl;
    // 	cout<<" N lines: "<<n.nLines<<endl;
    // 	cout<<" Points: "<<n.iPoint<<endl;
    // 	cout<<" Free legs when assigning: "<<n.nFreeLegsWhenAssigning<<endl;
    // 	cout<<" N poss: "<<n.nPoss<<endl;
    //   }
    
    possTable.resize(2*nnAss.size());
    
    /// Tensor product of all assignment heads and tail case
    vector<int64_t> curr(2*nnAss.size());
    
    /// Last assignemnt before overflow
    vector<int64_t> last(2*nnAss.size());
    
    for(int iNnAss=0;iNnAss<(int)nnAss.size();iNnAss++)
      {
	/// Nonnull ass
	const NnAss<S>& a=
	  nnAss[iNnAss];
	
	/// Index of the assignment
	const int& iAss=
	  a.iAss;
	
	/// Number of legs to assign is equal to the number of lines of the assignment
	const S& nLegsToAss=
	  a.nLines;
	
	/// Number of free legs when assigning the head
	const S& nFreeLegsFrom=
	  nFreeLegsWhenAssigning[iAss][FROM];
	
	/// Last possible assignment on the head
	const int64_t lastPossFrom=
	  lastCombination(nLegsToAss,nFreeLegsFrom);
	
	/// Store the starting side combination
	for(int64_t possFrom=firstCombination(nLegsToAss,nFreeLegsFrom);
	    possFrom<=lastPossFrom;
	    possFrom=nextCombination(possFrom))
	  possTable[2*iNnAss+FROM].push_back(decryptCombination(nLegsToAss,nFreeLegsFrom,possFrom));
	
	/// Number of free legs when assigning the tail
	const S& nFreeLegsTo=
	  nFreeLegsWhenAssigning[iAss][TO];
	
	/// Last possible assignment on the tail
	const int64_t lastPossTo=
	  lastDisposition(nLegsToAss,nFreeLegsTo);
	
	/// Store the ending side disposition
	for(int64_t possTo=firstDisposition(nLegsToAss,nFreeLegsTo);
	    possTo<=lastPossTo;
	    possTo=nextDisposition(possTo))
	  possTable[2*iNnAss+TO].push_back(decryptDisposition(nLegsToAss,nFreeLegsTo,possTo));
      }
    
    possibilitiesLooper=
      make_unique<Digits<S>>(fillVector<S>(possTable.size(),[this](const S& i)
							   {
							     return
							       possTable[i].size();
							   }));
  }
  
  WicksFinder(const vector<S>& nLegsPerPoint,const Assignment<S>& ass) :
    nLegsPerPoint(nLegsPerPoint),
    nPoints(nLegsPerPoint.size()),
    nLegs(summatorial(nLegsPerPoint)),
    nLines(nLegs/2),
    pointOfLeg(fillVector<S>(nLegs,[&](const S ileg)
			       {
				 /// Index of the output point
				 S iPoint=
				   0;
				 
				 /// Sum of all legs of points [0,ipoint)
				 S sumLegs=
				   0;
				 
				 while(ileg>=sumLegs+nLegsPerPoint[iPoint])
				   sumLegs+=nLegsPerPoint[iPoint++];
				 
				 return
				   iPoint;
			       })),
    ass(ass),
    nLegsPermPerPoint(transformVector(nLegsPerPoint,factorial<int64_t>)),
    nPermPerAss(transformVector(ass,factorial<int64_t>)),
    nLegsPermAllPoints(productorial(nLegsPermPerPoint)),
    nFreeLegsWhenAssigning(getNFreeLegsWhenAssigning()),
    nPermAllAss(productorial(nPermPerAss)),
    nLegsBefPoint(fillVector<S>(nLegsPerPoint.size(),[&](const S& iPoint)
				  {
				    /// Result
				    S res=
				      0;
				    
				    for(S jPoint=0;jPoint<iPoint;jPoint++)
				      res+=nLegsPerPoint[jPoint];
				    
				    return
				      res;
				  }))
  {
    reset();
    
    // cout<<" ANNA propStr: "<<nLegsPerPoint<<endl;
    // cout<<" ANNA ass: "<<ass<<endl;
    // cout<<" ANNA nLegsPermPerPoint: "<<nLegsPermPerPoint<<endl;
    // cout<<" ANNA nFreeLegsWhenAssigning: "<<nFreeLegsWhenAssigning<<endl;
  }
};

/// Compute the number of all Wick contractions
template <typename S>
int64_t computeNTotWicks(const vector<Assignment<S>>& allAss,const vector<S>& nPoints,const bool verbose=true)
{
  /// Result returned
  int64_t nTotWicks=
    0;
  
  if(verbose)
    COUT<<"List of all assignments: "<<endl;
  
  for(auto& ass : allAss)
    {
      const int64_t nWicks=
	WicksFinder<S>(nPoints,ass).nAllWickContrs(false);
      
      if(verbose)
	COUT<<" "<<ass<<" nWick: "<<nWicks<<endl;
      
      nTotWicks+=
	nWicks;
    }
  
  return
    nTotWicks;
}

#endif
