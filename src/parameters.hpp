/*
 * parameter.hpp
 *
 *  Created on: May 9, 2015
 *      Author: dmarce1
 */

#ifndef PARAMETER_HPP_
#define PARAMETER_HPP_

#include "RKDGFMM.hpp"
#include <string>

class parameters {
private:
	std::string filename;
	bool find_parameter( std::string& rstr, const std::string&) const;
public:
	parameters(const char*);
	bool get_parameter( real&, const char* ) const;
	bool get_parameter( integer&, const char* ) const;
	bool get_parameter( std::string&, const char* ) const;
};

#endif /* PARAMETER_HPP_ */
