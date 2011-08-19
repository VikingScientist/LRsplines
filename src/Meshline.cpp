
#include "Meshline.h"

namespace LR {

Meshline::Meshline(bool u_line, double const_par, double start, double stop, int multiplicity) {
	u_line_       =  u_line        ;
	const_par_    =  const_par     ;
	start_        =  start         ;
	stop_         =  stop          ;
	multiplicity_ =  multiplicity  ;
}


} // end namespace LR


