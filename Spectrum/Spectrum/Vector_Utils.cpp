#ifndef ATTACH_H
#include "Attach.h"
#endif

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