#ifndef _TOOLS_HPP
#define _TOOLS_HPP

#include <functional>
#include <vector>

using namespace std;

/// Compute the square
template <typename T>
T sqr(const T& t)
{
  return t*t;
}

/// Factorial
template <typename T>
T factorial(const T& n)
{
  /// Result
  T res=
    1;
  
  for(T i=2;i<=n;i++)
    res*=
      i;
  
  return res;
}

/// Transforms a vector using a function f
template <typename Tin,
	  typename F>
auto transformVector(const vector<Tin>& in,F f)
{
  /// Result
  vector<decltype(f(in[0]))> out;
  
  out.resize(in.size());
  
  for(int i=0;i<(int)in.size();i++)
    out[i]=f(in[i]);
  
  return out;
}

/// Fills a vector using a function f
template <typename T,
	  typename F>
vector<T> fillVector(const int n,F f)
{
  /// Result
  vector<T> out(n);
  
  for(int i=0;i<n;i++)
    out[i]=
      f(i);
  
  return out;
}

/// Transforms a vector using a function f
template <typename Tin,
	  typename F,
	  typename Tout>
Tout reduceVector(const vector<Tin>& in,F f,Tout res)
{
  for(int i=0;i<(int)in.size();i++)
    res=f(res,in[i]);
  
  return res;
}

/// Product of all elements of a vector
template <typename T>
T productorial(const vector<T>& in)
{
  return
    reduceVector(in,multiplies<int64_t>(),(T)1);
}

/// Sum of all elements of a vector
template <typename T>
T summatorial(const vector<T>& in)
{
  return
    reduceVector(in,plus<int64_t>(),(T)0);
}

/// Decompose a number in its digit following a multi-basis
template <typename T>
vector<T> decomposeNumber(T num,const vector<T>& basis)
{
  /// Result
  vector<T> out(basis.size());
  
  for(int i=basis.size()-1;i>=0;i--)
    {
      out[i]=
	num%basis[i];
      
      num/=
	basis[i];
    }
  
  return out;
}


/// Decompose a number in its digit following a fixed basis n
template <typename T>
vector<T> decomposeNumber(T num,const int nDigits,const int base)
{
  return decomposeNumber(num, fillVector<T>(nDigits,[base](const int i){return base;}));
}

/// Takes a permutation id and returns its assigned choice
template <typename T>
vector<int> decryptPermutation(const int n,T in)
{
  /// Returned permutation element
  vector<int> out(n,-1);
  
  /// Mask to determine the digit in the permutation
  int mask=
    n;
  
  for(int i=0;i<n;i++)
    {
      // Choice of the i-th element
      int c=
	in%mask;
      
      /// Poisition is initially 0
      int p=
	0;
      
      // Find the c-th non null and not used
      while(c!=0 or out[p]>=0)
	{
	  // If current is free, decrease number to move for the choice
	  if(out[p]<0)
	    c--;
	  
	  // Increment target position
	  p++;
	}
      
      // Mark down
      out[p]=i;
      
      in/=mask;
      
      mask--;
    }
  
  return out;
}

#endif
