#ifndef _LR_STREAMABLE_H
#define _LR_STREAMABLE_H

#include <iostream>

namespace LR {

class Streamable {
public:
	Streamable() {};
	virtual ~Streamable() {};

	virtual void read(std::istream& is)        = 0;
	virtual void write(std::ostream& os) const = 0;

};

inline std::istream& operator>>(std::istream& is, LR::Streamable& obj) {
	obj.read(is);
	return is;
}

inline std::ostream& operator<<(std::ostream& os, const LR::Streamable& obj) {
	obj.write(os);
	return os;
}

} // namespace LR

#endif
