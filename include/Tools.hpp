#ifndef _TOOLS_HPP
#define _TOOLS_HPP

#include <bitset>
#include <chrono>
#include <functional>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include <mpi.h>

using namespace std;

#define RED "\x1b[31m"
#define GREEN "\x1b[32m"
#define DEFAULT "\x1b[39m"

/// Compute the square
template <typename T>
T sqr(const T& t)
{
  return t*t;
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
template <typename Tin,
	  typename Tout>
vector<Tout> decomposeNumber(Tin num,const vector<Tout>& basis)
{
  /// Result
  vector<Tout> out(basis.size());
  
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
template <typename Tin,
	  typename Tout>
vector<Tout> decomposeNumber(Tin num,const int nDigits,const Tout base)
{
  return decomposeNumber(num,fillVector<Tout>(nDigits,[base](const Tout i){return base;}));
}

/// Prints an array
template <typename T,
	  uint64_t N>
ostream& operator<<(ostream& os,const array<T,N>& a);

/// Prints a vector
template <typename T>
ostream& operator<<(ostream& os,const vector<T>& a);

/// Prints using automatic range
template <typename T>
ostream& rangePrint(ostream& os,const T& a)
{
  // Open the bracket
  os<<"(";
  
  // Header
  bool h=
    true;
  
  for(auto& i : a)
    {
      if(h)
	h=false;
      else
	os<<",";
      
      os<<i;
    }
  
  os<<")";
  
  return
    os;
}

/// Prints an array
template <typename T,
	  uint64_t N>
ostream& operator<<(ostream& os,const array<T,N>& a)
{
  return
    rangePrint(os,a);
}

/// Prints a vector
template <typename T>
ostream& operator<<(ostream& os,const vector<T>& a)
{
  return
    rangePrint(os,a);
}

/// Cast to bitset
template <typename T,
	  int N=8*sizeof(T)>
bitset<N> bitRep(const T& t)
{
  /// Bitset
  return
    *(reinterpret_cast<const bitset<N>*>(&t));
}

/// Machine clock type
using Clock=
		   std::chrono::steady_clock;

/// Instant type, defining a moment in time
using Instant=
		   std::chrono::time_point<Clock>;

/// Difference of time between two moments
using Duration=
		   decltype(Instant{}-Instant{});

/// Get current time
inline Instant takeTime()
{
  return
    Clock::now();
}

/// Convert duration into seconds
template <typename O=double>                    // Output type
double durationInSec(const Duration& duration) ///< Input duration
{
  return
    std::chrono::duration<O>(duration).count();
}

template <typename K,typename V>
map<K,V> allReduceMap(const map<K,V>& in)
{
  map<K,V> out;
  
  int nRanks;
  MPI_Comm_size(MPI_COMM_WORLD,&nRanks);
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  for(int jRank=0;jRank<nRanks;jRank++)
    {
      int size=
	in.size();
      
      MPI_Bcast(&size,1,MPI_INT,jRank,MPI_COMM_WORLD);
      for(auto p : in)
	{
	  MPI_Bcast(&p,sizeof(p),MPI_CHAR,jRank,MPI_COMM_WORLD);
	  out[p.first]+=p.second;
	}
    }
  
  return
    out;
}

#endif
