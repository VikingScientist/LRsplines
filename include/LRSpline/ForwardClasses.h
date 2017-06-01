#ifndef FORWARDCLASSES_H
#define FORWARDCLASSES_H

#ifdef HAS_BOOST
	#include <boost/shared_ptr.hpp>
#else
#include <memory>
#endif

namespace LR {
	class Element;
	class Meshline;
	class MeshRectangle;
	class Basisfunction;
	class LRSplineSurface;
	class LRSplineVolume;
	class LRSpline;
}

#ifdef HAS_BOOST
	typedef boost::shared_ptr<LR::Element>       ElementPointer;
	typedef boost::shared_ptr<LR::Basisfunction> BasisPointer;
	typedef boost::shared_ptr<LR::MeshRectangle> MeshRectPointer;
	typedef boost::shared_ptr<LR::Meshline>      MeshlinePointer;
#else
	typedef std::shared_ptr<LR::Element>         ElementPointer;
	typedef std::shared_ptr<LR::Basisfunction>   BasisPointer;
	typedef std::shared_ptr<LR::MeshRectangle>   MeshRectPointer;
	typedef std::shared_ptr<LR::Meshline>        MeshlinePointer;
#endif


#endif
