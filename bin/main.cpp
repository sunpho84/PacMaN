#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include "Assignment.hpp"

/// Compute the square
template <typename T>
T sqr(const T& t)
{
  return t*t;
}

int main(int narg,char **arg)
{
  /// Defines the N-Point function
  const vector<int> nPoint=
    {5,5,2,4,4,2};
  
  /// Finder of all assignments
  AssignmentsFinder assignmentsFinder(nPoint);
  
  /// All assigments
  vector<Assignment> allAss=
    assignmentsFinder.getAllAssignements();
  
  drawAllAssignments(cerr,allAss,nPoint);
  
  return 0;
}
