#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H

// Declaration of a name space in which several useful std::vector related functions are defined
// R. Sheehan 7 - 12 - 2018

namespace vecut {
	void read_into_vector(std::string &filename, std::vector<double> &data, int &n_pts, bool loud = false);
	void write_into_file(std::string &filename, std::vector<double> &data, bool loud = false); 

	void wrap_around_conversion(std::vector<double> &data, bool ToStandard = true); 

	void read_into_matrix(std::string& filename, std::vector<std::vector<double>>& data, int& n_rows, int& n_cols, bool loud = false);
}

#endif
