#include <iostream>
#include <fstream>
#include "LRSpline/LRSplineSurface.h"

using namespace std;
using namespace LR;

int main(int argc, char **argv) {
  // create objects
  LRSplineSurface lr1(5,5,4,4);
  LRSplineSurface lr2(5,5,4,4);

  // do some refinements
  lr1.refineElement(0);
  lr1.refineElement(4);
  lr1.refineElement(6);
  lr1.refineElement(4);
  vector<double> knots = lr1.getEdgeKnots(WEST, true);
  cout << "Edge nodes on WEST side:\n";
  for(auto d : knots)
    cout << d << endl;
  // manually reverse knots (0,1) -> (1,0)
  for(int i=0; i<knots.size(); i++)
    knots[i] = 1-knots[i];

  // match on other side
  lr2.matchParametricEdge(NORTH, knots, true);
  cout << lr2 << endl;

  // dump files for easy inspection
  ofstream out1("lr1.eps");
  lr1.writePostscriptMesh(out1);
  ofstream out2("lr2.eps");
  lr2.writePostscriptMesh(out2);
  out1.close();
  out2.close();
  
} 
