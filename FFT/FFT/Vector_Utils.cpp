#include "pch.h" // use stdafx.h in Visual Studio 2017 and earlier

void vecut::read_into_vector(std::string &filename, std::vector<double> &data, int &n_pts, bool loud)
{
	// read a single column of data from a file into a vector
	// It is assumed that data is empty when the values are being read in
	// R. Sheehan 11 - 9 - 2017

	try {
		std::ifstream the_file;
		the_file.open(filename, std::ios_base::in);

		if (the_file.is_open()) {

			if (loud) std::cout << filename << " opened for reading\n";

			if (!data.empty()) data.clear(); // make sure that data is empty when values are being read in

			double value;
			n_pts = 0;
			while (the_file >> value) {
				data.push_back(value);
				n_pts++;
			}

			if (loud) std::cout << template_funcs::toString(n_pts) << " data were read from " << filename << "\n";

			the_file.close();
		}
		else {
			std::string reason;
			reason = "Error: void read_into_vector(std::string &filename, std::vector<double> &data)\n";
			reason += "Cannot open: " + filename + "\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void vecut::write_into_file(std::string &filename, std::vector<double> &data, bool loud)
{
	// write a single column of data to a new file

	try {

		if ( !data.empty() && filename != empty_str && useful_funcs::valid_filename_length(filename) ) {

			std::ofstream write(filename, std::ios_base::out, std::ios_base::trunc);

			if (write.is_open()) {

				for (size_t k = 0; k < data.size(); k++) {
					write << std::setprecision(10) << data[k] << "\n"; 
				}
			
				write.close(); 
			}
			else {
				std::string reason; 
				reason = "Error: void vecut::write_into_file(std::string &filename, std::vector<double> &data, bool loud)\n";
				reason += "Could not open file: " + filename + "\n"; 
				throw std::runtime_error(reason); 
			}
		}
		else {
			std::string reason;
			reason = "Error: void vecut::write_into_file(std::string &filename, std::vector<double> &data, bool loud)\n";
			reason += "Filename: " + filename + " is not valid\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error &e) {
		std::cerr << e.what(); 
	}
}

void vecut::wrap_around_conversion(std::vector<double>& data, bool ToStandard)
{
	// exchange the ranges data[_First, _Mid) and data[_Mid, _Last) in place
	// data is a vector of non-zero length
	// ToStandard == True implies data is in Wrap-Around order and you want to convert to Standard Ordering
	// ToStandard == False implies data is in Standard order and you want to convert to Wrap-Around Ordering
	// 
	// convert the elements in the vector data from wrap-around order to standard order
	// take elements from positions [size()/2, size()-1] and move them to positions [0, -1 + size()/2]
	// take elements from positions [0, -1 + size()/2] and move them to positions [size()/2, size()-1]
	// 
	// it would be nice if you could do this in place and it turns out that you can
	// While looking for a generic function to convert the output of the FFT algorithm from wrap-around to standard ordering
	// I encountered the std::rotate method see https://cplusplus.com/reference/algorithm/rotate/ for details
	// rotate the order of the elements in the range [first,last), in such a way that the element pointed by middle becomes the new first element
	// in other words exchange the ranges [_First, _Mid) and [_Mid, _Last)
	// 
	// Odd-ness or Even-ness of data.size() must be accounted for
	// 
	// R. Sheehan 22 - 8 - 2022

	try {
		if (!data.empty()) {

			size_t length = data.size();
			size_t half_length = length / 2;
			size_t mid_point = ToStandard ? (length % 2 == 0 ? half_length : half_length + 1) : half_length;

			std::rotate(data.begin(), data.begin() + mid_point, data.end()); // will convert between wrap-around and standard order in place
		}
		else {
			std::string reason;
			reason = "Error: void vecut::wrap_around_conversion(std::vector<double>& data)\n";
			reason += "data is empty\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void vecut::read_into_matrix(std::string& filename, std::vector<std::vector<double>>& data, int& n_rows, int& n_cols, bool loud)
{
	// read an array of data from a file
	// store the data in a matrix of size n_rows * n_cols
	// R. Sheehan 18 - 12 - 2018

	try {
		std::ifstream the_file;
		the_file.open(filename, std::ios_base::in);

		if (the_file.is_open()) {

			if (loud) std::cout << filename << " opened for reading\n";

			if (!data.empty()) data.clear(); // make sure that data is empty when values are being read in

			bool has_header = false;
			std::string line, item;

			char endline = '\n';
			char tab_token = '\t';
			char comma_token = ',';
			char sep_token;

			// Initialise the values to zero
			n_rows = 0; n_cols = 0;

			// Determine which token separates data in the file
			the_file.seekg(1, std::ios::beg); // move to the start of the file
			std::getline(the_file, line, endline);
			if (line.find(tab_token) != std::string::npos) sep_token = tab_token;
			if (line.find(comma_token) != std::string::npos) sep_token = comma_token;

			if (loud) std::cout << filename << " uses " << sep_token << " as a token\n";

			// Determine if first line is a file header
			if (isalpha(line[2])) has_header = true;

			// Count the number of rows and columns
			// This only seems to work when the data are separated by ',' also works for '\t' and ' '
			// http://www.cplusplus.com/reference/string/string/getline/
			// getline (istream& is, string& str, char delim)
			// Extracts characters from is and stores them into str until the delimitation character delim is found
			the_file.seekg(has_header ? 1 : 0, std::ios::beg); // move to the start of the file
			while (std::getline(the_file, line, endline)) {
				n_rows++;
				std::istringstream linestream(line);
				if (n_rows == 1) {
					while (std::getline(linestream, item, sep_token)) {
						n_cols++;
					}
				}
			}

			if (loud) std::cout << filename << " contains " << n_rows << " rows and " << n_cols << " columns\n";

			if (n_rows > 1 && n_cols > 0) {
				// Allocate the memory required to store the data
				data.resize(n_rows);
				for (size_t k = 0; k < data.size(); k++) {
					data[k].resize(n_cols, 0.0);
				}

				the_file.clear(); // empty a buffer, needed to ensure data can be read from the file
				the_file.seekg(has_header ? 1 : 0, std::ios::beg); // move to the start of the file

				int i, j;

				i = 0;
				while (std::getline(the_file, line, endline)) {
					std::istringstream linestream(line);
					j = 0;
					while (std::getline(linestream, item, sep_token)) {
						data[i][j] = atof(item.c_str());
						j++;
					}
					i++;
				}

				the_file.close();
			}
			else {
				std::string reason;
				reason = "Error: void read_into_matrix(std::string &filename, std::vector<std::vector<double>> &data, int &n_rows, int &n_cols, bool loud = false)\n";
				reason = filename + " contains no data\n";
				reason += "n_rows: " + template_funcs::toString(n_rows) + ", n_cols: " + template_funcs::toString(n_cols) + "\n";
				throw std::runtime_error(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void read_into_matrix(std::string &filename, std::vector<std::vector<double>> &data, int &n_rows, int &n_cols, bool loud = false)\n";
			reason += "Cannot open: " + filename + "\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
	catch (std::runtime_error& e) {
		std::cerr << e.what();
	}
}