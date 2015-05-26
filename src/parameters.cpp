/*
 * parameter.cpp
 *
 *  Created on: May 9, 2015
 *      Author: dmarce1
 */

#define LINE_LEN 256
#include "parameters.hpp"
#include <cstring>

bool parameters::find_parameter(std::string& rstr, const std::string& tmp_parameter) const {
	bool rc = false;
	const auto plen = tmp_parameter.size();
	std::string parameter = tmp_parameter;
	std::transform(parameter.begin(), parameter.end(), parameter.begin(), tolower);

	FILE* fp = fopen(filename.c_str(), "rt");
	char _line[LINE_LEN];
	std::string line;
	if (fp == nullptr) {
		printf("Unopen to open parameter file %s, aborting...\n", filename.c_str());
		abort();
	}
	while (1) {
		fgets(_line, LINE_LEN, fp);
		line = _line;
		std::transform(line.begin(), line.end(), line.begin(), tolower);
		const char* ptr = line.c_str();
		while (isspace(*ptr)) {
			++ptr;
		}
		if (strncmp(ptr, parameter.c_str(), plen) == 0) {
			rstr = ptr + plen;
			rc = true;
		}
	}
	fclose(fp);
	return rc;
}

parameters::parameters(const char* str) {
	filename = str;
}
bool parameters::get_parameter(real& f, const char* name) const {
	bool rc;
	std::string value;
	rc = find_parameter(value, name);
	if (rc) {
		f = std::atof(value.c_str());
	}
	return rc;
}

bool parameters::get_parameter(integer& i, const char* name) const {
	bool rc;
	std::string value;
	rc = find_parameter(value, name);
	if (rc) {
		i = std::atoi(value.c_str());
	}
	return rc;
}

bool parameters::get_parameter(std::string& s, const char* name) const {
	bool rc;
	std::string value;
	rc = find_parameter(value, name);
	if (rc) {
		s = value;
	}
	return rc;
}
