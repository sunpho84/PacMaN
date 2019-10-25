#ifndef _COMBINATORIAL_HPP
#define _COMBINATORIAL_HPP

#include <cstdint>
#include <vector>

using namespace std;

/// Partition of a number
template <typename S>
using Partition=
  vector<S>;

/// List all partitioning of the number m
template <typename S>
vector<Partition<S>> listAllPartitioningOf(const S& m)
{
  /// List to be returned
  vector<Partition<S>> out;
  
  /// Position
  S k=
    0;
  
  /// Current vector
  vector<S> p(m,0);
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
      S remVal=
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

/// Factorial
template <typename T=int64_t>
T factorial(const int& n)
{
  /// Result
  T res=
    1;
  
  for(T i=2;i<=n;i++)
    res*=
      i;
  
  return res;
}

/// Ratio of factorials num!/den!
template <typename T=int64_t>
T factorialsRatio(const int& num,const int& den)
{
  /// Result
  T res=
    1;
  
  for(T i=std::max(2,den+1);i<=num;i++)
    res*=
      i;
  
  return res;
}

/// Netwon binomial (n m)
template <typename T=int64_t>
T newtonBinomial(const int& n,const int& m)
{
  if(n<m)
    return 0;
  else
    return
      factorialsRatio<T>(n,max(m,n-m))/factorial<T>(min(m,n-m));
}

/// Returns the number of permutations of n
template <typename T=int64_t>
T nPermutations(const int& n)
{
  return
    factorial<T>(n);
}

/// Takes the id disposition of nObj numbers into nSlots slots, and returns its assigned choice
template <typename T,
	  typename S>
vector<S> decryptDisposition(const S& nObj,const int& nSlots,T iDisp)
{
  /// Choice at each turn
  vector<S> choice(nObj);
  
  /// Mask to determine the digit in the permutation
  S mask=
    nSlots-nObj+1;
  
  for(S iObj=nObj-1;iObj>=0;iObj--)
    {
      // Choice of the i-th object
      choice[iObj]=
	iDisp%mask;
      
      iDisp/=mask;
      
      mask++;
    }
  
  // Increment chosen one
  for(S iObj=nObj-2;iObj>=0;iObj--)
    for(S jObj=iObj+1;jObj<nObj;jObj++)
      if(choice[jObj]>=choice[iObj])
	choice[jObj]++;
  
  return choice;
}

/// Takes the id disposition of n numbers into n slots,and returns its assigned choice
template <typename T,
	  typename S>
vector<S> decryptPermutation(const int& n,const T& in)
{
  return
    decryptDisposition<T,S>(n,n,in);
}

/// Returns the number of dispositions (ordered choices) of nObj objects into nSlots slots
template <typename T=int64_t,
	  typename S>
T nDispositions(const S& nObj,const S& nSlots)
{
  return
    factorialsRatio<T>(nSlots,nSlots-nObj);
}

/// Returns the first dispositions of nObj objects into nSlots slots
template <typename T=int64_t,
	  typename S>
T firstDisposition(const S& nObj,const S& nSlots)
{
  return
    0;
}

/// Returns the next disposition
template <typename T>
T nextDisposition(const T& x)
{
  return
    x+1;
}

/// Returns the last dispositions of nObj objects into nSlots slots
template <typename T=int64_t,
	  typename S>
T lastDisposition(const S& nObj,const S& nSlots)
{
  return
    nDispositions<T,S>(nObj,nSlots)-1;
}

/// Returns the number of combinations (unordered choices) of nObj objects into nSlots slots
template <typename T=int64_t,
	  typename S>
T nCombinations(const S& nObj,const S& nSlots)
{
  return
    newtonBinomial<T>(nSlots,nObj);
}

/// Takes the id combination of nObj numbers into nSlots slots, and returns its assigned choice
template <typename T,
	  typename S>
vector<S> decryptCombination(const S& nObj,const S& nSlots,T iCombo)
{
  /// Result
  vector<S> out(nObj);
  
  /// Index of object
  S iObj=
    0;
  
  /// Index of the slot
  S iSlot=
    0;
  
  while(iObj<nObj and iSlot<nSlots)
    {
      if((iCombo>>iSlot)&1)
	out[iObj++]=
	  iSlot;
      iSlot++;
    }
  
  return out;
}

// Find next k-combination
template <typename T>
T nextCombination(const T& extX)
{
  /// Unsigned version of T
  using U=
    make_unsigned_t<T>;
  
  /// Local unsigned copy
  const U x=
    static_cast<U>(extX);
  
  /// Extract rightmost bit 1
  const U u=
    x & -x;
  
  /// Set last non-trailing bit 0, and clear to the right
  const U v=
    u+x;
  
  return
    static_cast<T>(v+(((v^x)/u)>>2));
}

/// First combination of nObj objects
template <typename T=int64_t,
	  typename S>
T firstCombination(const S& nObj,const S& nSlots=0)
{
  return
    (static_cast<T>(1)<<nObj)-1;
}

/// Last combination of nObj objects into nSlots slots
template <typename T=int64_t,
	  typename S>
T lastCombination(const S& nObj,const S& nSlots)
{
  return
    ((static_cast<T>(1)<<nSlots)-1)^
    ((static_cast<T>(1)<<(nSlots-nObj))-1);
}

/// Loops on all combinations of nSlots
template <typename F,
	  typename T=int64_t>
void forAllCombinations(const int& nObj,const int& nSlots,F f)
{
  if(nObj>0 and nObj<=nSlots)
    {
      /// Last combination
      T last=
	lastCombination(nObj,nSlots);
      
      /// Running combo
      T combo=
	firstCombination(nObj);
      
      do
	{
	  f(combo);
	  
	  combo=
	    nextCombination(combo);
	}
      while(combo<=last);
    }
}

/////////////////////////////////////////////////////////////////

/// Represents a mixed basis number
template <typename S>
class Digits
{
public:
  
  /// Base of each digit
  const vector<S> base;
  
  /// Number of digits
  S nDigits() const
  {
    return
      base.size();
  };
  
  /// Last digit
  const S lastDigit;
  
  /// Digits rerpresenting the number
  vector<S> digits;
    
  Digits(const vector<S>& base) :
    base(base),
    lastDigit(nDigits()-1),
    digits(nDigits())
  {
  }
  
  /// Set to a given number
  template <typename T>
  void setTo(T t)
  {
    for(S iDigit=nDigits()-1;iDigit>=0;iDigit--)
      {
	// Choice of the i-th digitect
	digits[iDigit]=
	  t%base[iDigit];
	
	t/=
	  base[iDigit];
      }
  }
  
  /// Loop on all numbers
  template <typename F>
  void forAllNumbers(F f)
  {
    // Reset
    setTo(0);
    
    /// Index of the running digit
    S iDigit=
      lastDigit;
    
    // Loops on all numbers
    while(iDigit>=0)
      {
	// Exec the function
	f(digits);
	
	while(iDigit>=0 and (++digits[iDigit])>=base[iDigit])
	  digits[iDigit--]=
	    0;
	
	// If not run on all
	if(iDigit>=0)
	  iDigit=
	    lastDigit;
      }
  }
};

#endif
