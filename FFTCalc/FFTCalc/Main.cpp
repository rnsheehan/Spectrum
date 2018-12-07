#ifndef ATTACH_H
#include "Attach.h"
#endif // !ATTACH_H

// The aim of this project is to implement the code needed to compute the FFT of an arbitrary data set
// The data I have in mind is that obtained from a high speed oscilloscope
// Generally, the data will be of the form (time_value, voltage_value)
// The code for computing the FFT is taken from NRinC, this code adapts an existing code to run using std::vector<double>
// One of the tasks required is to determine the SNR of the FFT from the trace data
// R. Sheehan 5 - 12 - 2018

// This version of the code uses the developed FFT code and will run from the command line
// It is this version of the code that will run from inside LabVIEW
// R. Sheehan 7 - 12 - 2018

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

			std::string timefile = argv[1];
			std::string spctfile = argv[2];

			if (useful_funcs::file_exists(timefile) && useful_funcs::file_exists(spctfile)) {

				std::vector<double> timedata; int ntimes = 0;
				
				std::vector<double> spctdata; int nspct = 0;

				vecut::read_into_vector(timefile, timedata, ntimes, true);

				vecut::read_into_vector(spctfile, spctdata, nspct, true);

				double delta_t = timedata[1] - timedata[0];

				unsigned long nn = nspct;

				fft calc; // this is the FFT calculation object

				clock_t start = clock(); 

				calc._four1(spctdata, nn); // compute the FFT of spctdata

				clock_t finish = clock(); 

				calc.output_data(spctdata, delta_t, spctfile); // output the absolute value of the FFT as well as the frequency space data

				std::cout << "FFT Calculation Complete\n"; 
				std::cout << "Time taken: " << static_cast<double>((finish - start) / CLOCKS_PER_SEC) << " seconds\n";
			}
			else {
				std::string reason;
				std::string prog_name = argv[0];
				reason = "Error: " + prog_name + ".exe\n";
				reason += "One or both of the files does not exist\n";
				if (!useful_funcs::file_exists(timefile)) reason += "File: " + timefile + " does not exist\n";
				if (!useful_funcs::file_exists(spctfile)) reason += "File: " + spctfile + " does not exist\n";
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
