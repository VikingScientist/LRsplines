#ifndef HASHSET_H
#define HASHSET_H

#include <list>
#include <map>
#include <cstddef>

/*!
	\brief HashSet iterator which allows for iteration over the HashSet class
	\details Internally, the iterator loops over two stl containers, the hashcodes which are stored as buckets in an stl::map
	         and secondly, an stl::list for hash codes which collide to the same bucket, but are distinct (does not pass equality-test)
*/

template<typename T>
class HashSet_iterator
         :public std::iterator<std::forward_iterator_tag,     // type of iterator
                               T,ptrdiff_t,T*,T&>             // Info about iterator
{

	typedef typename std::map<long, std::list<T> >::iterator       iter;
	typedef typename std::map<long, std::list<T> >::const_iterator citer;
	typedef typename std::list<T>::iterator                        list_iter;

public:

	//! \brief Default constructor
	HashSet_iterator() {
	}

	//! \brief Default constructor
	//! \param majorIter iterator position in the hashcode map
	//! \param subIter   iterator position in the linked list for non-unique hash codes
	//! \param majorEnd  iterator position to the end of the hashcode map
	HashSet_iterator(iter majorIter, list_iter subIter, iter majorEnd) {
		this->majorIter = majorIter;
		this->subIter   = subIter;
		this->majorEnd  = majorEnd;
	}

	//! \brief Dereferencing the iterator returns an object of class <T>
	T& operator*() const {
		return *subIter;
	}

	//! \brief Dereferencing the iterator returns an object of class <T>
	T* operator->() const {
		return &(*subIter);
	}

	HashSet_iterator& operator=(const HashSet_iterator &other) {
		majorEnd  = other.majorEnd ;
		majorIter = other.majorIter;
		subIter   = other.subIter  ;
		return *this;
	}

	HashSet_iterator& operator++() {
		subIter++;
		if(subIter == majorIter->second.end()) {
			majorIter++;
			if(majorIter != majorEnd)
				subIter = majorIter->second.begin();
		}
		return *this;
	}
	/*
	HashSet_iterator operator++(int i) {
		myIter += i;
	}
	*/
	bool equal(HashSet_iterator const& rhs) const {
		return majorIter == rhs.majorIter && subIter == rhs.subIter;
	}

private:
	iter      majorIter;
	iter      majorEnd;
	list_iter subIter;

};

/*!
	\brief const version of the HashSet iterator
*/
template<typename T>
class HashSet_const_iterator
         :public std::iterator<std::forward_iterator_tag,     // type of iterator
                               T,ptrdiff_t,const T*,const T&> // Info about iterator
{

	typedef typename std::map<long, std::list<T> >::const_iterator iter;
	typedef typename std::list<T>::const_iterator                  list_iter;

public:

	//! \brief Default constructor
	HashSet_const_iterator() {
	}

	//! \brief Default constructor
	//! \param majorIter iterator position in the hashcode map
	//! \param subIter   iterator position in the linked list for non-unique hash codes
	//! \param majorEnd  iterator position to the end of the hashcode map
	HashSet_const_iterator(iter majorIter, list_iter subIter, iter majorEnd) {
		this->majorIter = majorIter;
		this->subIter   = subIter;
		this->majorEnd  = majorEnd;
	}

	//! \brief Dereferencing the iterator returns an object of class <T>
	const T& operator*() const {
		return *subIter;
	}

	//! \brief Dereferencing the iterator returns an object of class <T>
	const T* operator->() const {
		return &(*subIter);
	}

	HashSet_const_iterator& operator=(const HashSet_const_iterator &other) {
		majorEnd  = other.majorEnd ;
		majorIter = other.majorIter;
		subIter   = other.subIter  ;
		return *this;
	}

	HashSet_const_iterator& operator++() {
		subIter++;
		if(subIter == majorIter->second.end()) {
			majorIter++;
			if(majorIter != majorEnd)
				subIter = majorIter->second.begin();
		}
		return *this;
	}
	/*
	HashSet_iterator operator++(int i) {
		myIter += i;
	}
	*/
	bool equal(HashSet_const_iterator const& rhs) const {
		return majorIter == rhs.majorIter && subIter == rhs.subIter;
	}

private:
	iter      majorIter;
	iter      majorEnd;
	list_iter subIter;

};

template<typename T>
inline bool operator!=(HashSet_iterator<T> const& lhs, HashSet_iterator<T> const& rhs)
{
	return !lhs.equal(rhs);
}

template<typename T>
inline bool operator==(HashSet_iterator<T> const& lhs, HashSet_iterator<T> const& rhs)
{
	return lhs.equal(rhs);
}

template<typename T>
inline bool operator!=(HashSet_const_iterator<T> const& lhs, HashSet_const_iterator<T> const& rhs)
{
	return !lhs.equal(rhs);
}

template<typename T>
inline bool operator==(HashSet_const_iterator<T> const& lhs, HashSet_const_iterator<T> const& rhs)
{
	return lhs.equal(rhs);
}


/*!
	\brief HashSet container which allows for quick sorting on a non-unique hashfunction, and only contains truly unique elements
	\details The container requires the class to implement the hashCode function which is used for sorting the results and allow
	         for access in logarithmic time. Where hash codes coincide, multiple elements are tested by a potentially more time consuming
	         equals-functions, and unique elements are stored in a linked list at this bucket. It is required that elements that pass
	         equality-test always produce the same hash code.
*/

template <class T>
class HashSet {

typedef typename std::map<long, std::list<T> >::iterator       iter;
typedef typename std::map<long, std::list<T> >::const_iterator citer;
typedef typename std::list<T>::iterator                        list_iter;

public:

	//! \brief Default constructor
	//! \details Creates an empty container

	HashSet() {
		numb      = 0;
	}

	//! \brief Default copy constructor
	//! \details Creates a copy of the HashSet.
	HashSet(const HashSet<T> &other) {
		data = other.data;
		numb = other.numb;
	}

	//! \brief insert an element in the container if it does not already exist
	//! \param obj the element to insert
	//! \details on the occasion that the element do already exist in the set, no action is taken.
	//!          Complexity: logarithmic in size, linear in hash collisions
	void insert(const T &obj) {
		long hc = obj->hashCode();
		iter it = data.find(hc);

		if(it == data.end()) {
			data[hc] = std::list<T>(1,obj);
			numb++;
		} else {
			for(list_iter lit = it->second.begin(); lit != it->second.end(); lit++)
				if((*lit)->equals(*obj))
					return;
			data[hc].push_back(obj);
			numb++;
		}
	}

	//! \brief erase an element in the container if it does exist
	//! \param obj the element to remove
	//! \returns 0 if no elements were removed, 1 if it did exist and was successfully removed
	//! \details on the occasion that the element do already exist in the set, no action is taken.
	//!          Complexity: logarithmic in size, linear in hash collisions
	int erase(const T &obj) {
		long hc = obj->hashCode();
		iter it = data.find(hc);
		if(it == data.end())
			return 0;

		for(list_iter lit = it->second.begin(); lit != it->second.end(); lit++) {
			if(obj->equals(**lit)) {
				it->second.erase(lit);
				if(it->second.size() == 0) {
					data.erase(it);
				}
				numb--;
				return 1;
			}
		}

		return 0;
	}

	//! \brief Searches the container for an element and returns an iterator to it if found, otherwise it returns an iterator to map::end (the element past the end of the container).
	//! \param obj The element to search for
	//! \details Complexity: logarithmic in size, linear in hash collisions
	HashSet_iterator<T> find(const T &obj) {
		long hc = obj->hashCode();
		iter it = data.find(hc);
		if(it == data.end())
			return end();
		for(list_iter lit = it->second.begin(); lit != it->second.end(); lit++)
			if(obj->equals(**lit))
				return HashSet_iterator<T>(it, lit, data.end());

		return end();
	}

	//! \brief returns the first element, and removes this from the container
	T pop() {
		if(numb == 0)
			return NULL;

		iter it = data.begin();
		T ans   = it->second.front();
		it->second.pop_front();
		if(it->second.size() == 0)
			data.erase(it);

		numb--;
		return ans;
	}

	//! \brief clears the container
	void clear() {
		for(iter it = data.begin(); it != data.end(); it++)
			it->second.clear();
		data.clear();
		numb = 0;
	}

	//! \brief returns the number of unique hash codes in this container
	//! \details this is the number of elements that produce distinct value upon calling T::hashFunction()
	int uniqueHashCodes() const {
		return data.size();
	}

	//! \brief returns the number of unique elements in the container
	//! \details this is the number of elements that produce false upon calling T::equals()
	int size() const {
		return numb;
	}

	//! \brief iterator to the beginning of the container.
		//! \details dereferencing the iterator returns an object of class <T>
	HashSet_const_iterator<T> begin() const {
		return HashSet_const_iterator<T>(data.begin(),
		                                 (numb==0)?dummyLast.end():data.begin()->second.begin(),
		                                 data.end());
	}

	//! \brief iterator to one past the last element
	//! \details dereferencing the iterator returns an object of class <T>
	HashSet_const_iterator<T> end() const {
                citer end = data.end();
                if (!data.empty())
                  end--;
		return HashSet_const_iterator<T>(data.end(),
		                                 (numb==0)?dummyLast.end():end->second.end(),
		                                 data.end());
	}


	//! \brief iterator to the beginning of the container.
	//! \details dereferencing the iterator returns an object of class <T>
	HashSet_iterator<T> begin() {
		return HashSet_iterator<T>(data.begin(),
		                           (numb==0)?dummyLast.end():data.begin()->second.begin(),
		                           data.end());
	}

	//! \brief iterator to one past the last element
	//! \details dereferencing the iterator returns an object of class <T>
	HashSet_iterator<T> end() {
                iter end = data.end();
                if (!data.empty())
                  end--;
		return HashSet_iterator<T>(data.end(),
		                           (numb==0)?dummyLast.end():end->second.end(),
		                           data.end());
	}

private:
	std::map<long, std::list<T> > data;
	std::list<T> dummyLast;
	int  numb;

};


#endif
