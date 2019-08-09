#ifndef USEFUL_H
#define USEFUL_H

// Library of functions that are very useful
// R. Sheehan 4 - 7 - 2011

namespace useful_funcs{
	
	std::string TheTime();

	void exit_failure_output(std::string reason);

	void remove_substring(std::string &the_string, std::string the_sub_string); 

	void create_directory(std::string &dir_name); 

	void set_directory(std::string &dir_name);

	void get_directory(); 

	bool valid_filename_length(const std::string &name); 

	unsigned long next_POT(double x); // This function converts a number x to the next highest power of two 

	bool is_POT(int x); // Test a number to see if it is a power of two
}

#endif