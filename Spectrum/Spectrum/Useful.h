#ifndef USEFUL_H
#define USEFUL_H

// Library of functions that are very useful
// R. Sheehan 4 - 7 - 2011

namespace useful_funcs{
	
	std::string TheTime();

	void exit_failure_output(std::string reason);

	void read_into_vector(std::string &filename, std::vector<double> &data, int &n_pts, bool loud = false); 

	unsigned long next_POT(double x); // This function converts a number x to the next highest power of two 

	bool is_POT(int x); // Test a number to see if it is a power of two
}

#endif