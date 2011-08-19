#ifndef ELEMENT_H
#define ELEMENT_H

namespace LR {

class Element {

public:
	Element(double start_u, double start_v, double stop_u, double stop_v);

private:
	double start_u_;
	double start_v_;
	double stop_u_;
	double stop_v_;

};

} // end namespace LR

#endif


