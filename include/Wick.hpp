#ifndef _WICK_HPP
#define _WICK_HPP

#include <memory>

#include "Assignment.hpp"

using namespace std;

/// A line connects two points
typedef array<int,2> Line;

/// Wick contraction is a list of lines
typedef vector<Line> Wick;

/// Creates all Wick contraction, given a n-point function and an assignment
class WicksFinder
{
  /// Number of legs in each point
  const vector<int> nLegsPerPoint;
  
  /// Number of points
  const int nPoints;
  
  /// Total number of legs
  const int nLegs;
  
  /// Total number of lines, equal to half the legs
  const int nLines;
  
  /// Point where the leg is attached
  const vector<int> pointOfLeg;
  
  /// Assignment to be considered
  const Assignment ass;
  
  /// Number of permutations of the legs of each point
  const vector<int64_t> nLegsPermPerPoint;
  
  /// Number of permutations of all lines of the given propagator
  /// assignments, for each point
  const vector<int64_t> nPermPerAss;
  
  /// Total number of permutations of all legs of all points
  const int64_t nLegsPermAllPoints;
  
  /// Number of free legs in the head and in the tail when assigning the (i,j) assignment
  const vector<array<int,2>> nFreeLegsWhenAssigning;
  
  /// Product of the number of all permutations of all lines of the
  /// given propagator assignment
  const int64_t nPermAllAss;
  
  /// Number of legs before the given point
  const vector<int> nLegsBefPoint;
  
  /// Compute the number of free legs on the head and tail of all
  /// assignment when assigning each of them in turn
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
    
    /// Considered assignment
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
  
  /// Precomputed list of all possible assignment
  vector<vector<vector<int>>> possTable;
  
  /// Non-null associations
  vector<NnAss> nnAss;
  
  /// Looper on all possibilities
  unique_ptr<Digits> possibilitiesLooper;
  
public:
  
  /// Return first Wick contraction
  Wick getFirst()
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
  
  /// Convert the digits of the Wick contraction id written in terms of digits into an actual Wick contraction
  Wick convertDigitsToWick(const vector<int>& wickDigits)
    const
  {
    /// Store wether the leg is assigned
    vector<int> legIsAss(nLegs,false);
    
    /// Store the assignment of the legs, in form of lines connecting two legs
    Wick lineAss(nLines);
    
    /// Index of the line to assign
    int iLineToAss=
      0;
    
    for(int iNnAss=0;iNnAss<(int)nnAss.size();iNnAss++)
      {
	/// Number of legs for this assignment
	const int nLegsPerAss=
	  nnAss[iNnAss].nLines;
	
	/// Legs connected by the lines of this assignment
	///
	/// We need to temporarily store this, because in the table
	/// of possibility they are counted as a whole block
	vector<array<int,2>> pointsLegAss(nLegsPerAss);
	
	for(int iLine=0;iLine<nLegsPerAss;iLine++)
	  for(int ft=0;ft<2;ft++)
	    {
	      // At first, set the leg to the number of legs preceeding
	      // the point (which is the lable of the first leg of the point)
	      int& l=
		pointsLegAss[iLine][ft]=
		nLegsBefPoint[nnAss[iNnAss].iPoint[ft]];
	      
	      // Then skip needed unassigned legs
	      int count=
		possTable[2*iNnAss+ft][wickDigits[2*iNnAss+ft]][iLine];
	      while(count>0 or legIsAss[l])
		{
		  if(not legIsAss[l])
		    count--;
		  l++;
		}
	    }
	
	// Now we mark all assigned
	for(int iLegPerAss=0;iLegPerAss<nLegsPerAss;iLegPerAss++)
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
    possibilitiesLooper->forAllNumbers([&,this](const vector<int>& wickDigits)
				       {
					 /// Line assigments
					 Wick lineAss=
					   convertDigitsToWick(wickDigits);
					 
					 f(lineAss);
				       });
  }
  
  /// Get the Wick contraction numberiWick
  Wick get(const int64_t& iWick)
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
	const NnAss& a=
	  nnAss[iNnAss];
	
	/// Index of the assignment
	const int& iAss=
	  a.iAss;
	
	/// Number of legs to assign is equal to the number of lines of the assignment
	const int& nLegsToAss=
	  a.nLines;
	
	/// Number of free legs when assigning the head
	const int& nFreeLegsFrom=
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
	const int& nFreeLegsTo=
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
      make_unique<Digits>(fillVector<int>(possTable.size(),[this](const int& i)
							   {
							     return possTable[i].size();
							   }));
  }
  
  WicksFinder(const vector<int>& nLegsPerPoint,const Assignment& ass) :
    nLegsPerPoint(nLegsPerPoint),
    nPoints(nLegsPerPoint.size()),
    nLegs(summatorial(nLegsPerPoint)),
    nLines(nLegs/2),
    pointOfLeg(fillVector<int>(nLegs,[&](const int ileg)
			       {
				 /// Index of the output point
				 int iPoint=
				   0;
				 
				 /// Sum of all legs of points [0,ipoint)
				 int sumLegs=
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
    nLegsBefPoint(fillVector<int>(nLegsPerPoint.size(),[&](const int iPoint)
				  {
				    /// Result
				    int res=
				      0;
				    
				    for(int jPoint=0;jPoint<iPoint;jPoint++)
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
int64_t computeNTotWicks(const vector<Assignment>& allAss,const vector<int>& nPoints,const bool verbose=true);

#endif
