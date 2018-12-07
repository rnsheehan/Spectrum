#ifndef ATTACH_H
#include "Attach.h"
#endif // !ATTACH_H

// The aim of this project is to implement the code needed to compute the FFT of an arbitrary data set
// The data I have in mind is that obtained from a high speed oscilloscope
// Generally, the data will be of the form (time_value, voltage_value)
// The code for computing the FFT is taken from NRinC, this code adapts an existing code to run using std::vector<double>
// One of the tasks required is to determine the SNR of the FFT from the trace data
// R. Sheehan 5 - 12 - 2018

int main(int argc, char *argv[])
{
	// run program from command line
	// assuming that inputs are going to be filename1 filename2
	// check if files exist
	// read data from both files
	// perform operation on data in files
	// output results of operation to two files
	try {
		if (argc > 1) {
			// List off the input parameters
			// Program needs 2 or more parameters to run, remember that the name of the program is also considered a parameter
			// argv[0] = program name

			std::cout << "Name of the program is " << argv[0] << ".exe\n";
			std::cout << argc - 1 << " parameters were input into the program\n";
			for (int count = 1; count < argc; count++) {
				std::cout << "argv[" << count << "] = " << argv[count] << "\n";
			}
			std::cout << "\n";

			// Run the calculation based on the input parameters
			// filename1 = argv[1]
			// filename2 = argv[2]

			std::string filename1 = argv[1];
			std::string filename2 = argv[2];

			if (useful_funcs::file_exists(filename1) && useful_funcs::file_exists(filename2)) {
				int size1 = 0, size2 = 0;
				std::vector<double> data1, data2;
				vecut::read_into_vector(filename1, data1, size1, true);
				vecut::read_into_vector(filename2, data2, size2, true);
			}
			else {
				std::string reason;
				std::string prog_name = argv[0];
				reason = "Error: " + prog_name + ".exe\n";
				reason += "One or both of the files does not exist\n";
				if (!useful_funcs::file_exists(filename1)) reason += "File: " + filename1 + " does not exist\n";
				if (!useful_funcs::file_exists(filename2)) reason += "File: " + filename2 + " does not exist\n";
				throw std::invalid_argument(reason);
			}
		}
		else {
			throw std::invalid_argument("Not enough arguments were input\n");
		}
	}
	catch (std::invalid_argument &e) {
		std::cerr << e.what();
	}

	return 0;
}
