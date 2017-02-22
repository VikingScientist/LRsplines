#include <iostream>
#include <fstream>
#include <vector>
#include <LRSpline/LRSplineSurface.h>

using namespace std;
using namespace LR;

int main(int argc, char **argv) {

 // cubic tensor spline with a 7x6 control-net
 // defaults to the identity-mapping:
 //    x(u,v) = u
 //    y(u,v) = v
	LRSplineSurface lr(7, 6, 4, 4);


  // do a surface evaluation
  double u,v;
	vector<double> xy_point(2);
  u = 0.77;
  v = 0.23;
	lr.point(xy_point, u,v);
	cout << "(u,v) = " <<           u << ", " <<           v <<  endl;
	cout << "(x,y) = " << xy_point[0] << ", " << xy_point[1] <<  endl;

  // print all basis functions (tensor product so far)
  for(Basisfunction* b : lr.getAllBasisfunctions())
    cout << *b << endl;

  // refine basis function #17
  lr.refineBasisFunction(17);

  // write lr mesh to an eps image-file
  ofstream outfile("mesh.eps");
  lr.writePostscriptMesh(outfile);
  outfile.close();

} 
