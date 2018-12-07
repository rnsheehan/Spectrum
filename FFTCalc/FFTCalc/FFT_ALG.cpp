#ifndef ATTACH_H
#include "Attach.h"
#endif

// Definitions for the fft class

fft::fft()
{
	// Default constructor
	arr_size = 0; 
}

// Methods
// Users can only access these methods because they have exception handling
void fft::_four1(std::vector<double> &data, unsigned long &nn, int isign, bool FORMAT_DATA)
{
	// four1 algorithm with exception handling
	// nn is the length of the original data set, it gets changed by the pad data algorithm if necessary

	// Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1
	// Replace data[1..2*nn] by nn times its inverse dicrete Fourier transform, if isign is input as -1
	// data is a complex array of length nn, or equivalently, a real array of length 2*nn
	// nn MUST be an integer power of 2

	// data must be in the form output by format_four1 before fft algorithm can be applied
	// padding of data can be applied after the fact

	try {
		
		if ( !data.empty() ) {

			if (FORMAT_DATA) format_four1(data);

			pad_data(data, nn, true); // data input into four1 must be in cmplx format

			if (useful_funcs::is_POT(nn) && abs(isign) == 1) {

				// nn must be a power of two
				// data must be in the form (real, imag, real, imag, real, imag, ....), with all imag = 0, hence length 2*nn
				// isign = 1 for FT
				four1(data, nn, isign);

				// scale by 1/N in the case of an inverse FT
				// Not yet clear what effect this has on other algorithms
				/*if(isign == IFT){
				for(int i=1; i<=data.n_elems(); i++){
				data[i] /= nn;
				}
				}*/
			}
			else {
				std::string reason;
				reason = "Error: void fft::pad_data(std::vector<double> &data, unsigned long &nn, bool CMPLX_ARR)\n";
				reason = "Input arrays do not have the correct dimension in fft::_four1\n";
				reason = reason + "nn = " + template_funcs::toString(nn) + "\n";
				reason = reason + "2*nn = " + template_funcs::toString(2 * nn) + "\n";
				if (useful_funcs::is_POT(nn)) {
					reason = reason + "nn is a power of 2\n";
				}
				else {
					reason = reason + "nn is not a power of 2\n";
				}
				reason = reason + "data.n_elems() = " + template_funcs::toString(data.size()) + "\n";
				reason = reason + "isign = " + template_funcs::toString(isign) + "\n";
				throw std::invalid_argument(reason);
			}
		}
		else {
			std::string reason;
			reason = "Error: void fft::pad_data(std::vector<double> &data, unsigned long &nn, bool CMPLX_ARR)\n";
			reason += "Input array is empty\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fft::output_data(std::vector<double> data, double pos_spac, std::string &filename, bool wrap_around)
{
	// output the computed FFT data to a file
	// How do you handle the position data if the size of the data has been expanded? 
	// Don't need it because you are working in frequency space. This means that you should output frequency data
	// for the human readable format. 

	// Need two options for output
	// Standard option is to output directly from the FFT algorithm, this outputs in wrap-around order
	// Also want an option to output data in a human readable / plottable format

	// this assumes that FFT of data has been computed

	try {
		if (pos_spac > 0.0 && !data.empty() && filename != empty_str && useful_funcs::valid_filename_length(filename)) {
		
			// remove .txt from end of filename
			useful_funcs::remove_substring(filename, dottxt);

			if (wrap_around) {
				// output data in standard wrap-around order
				// never really want to use this order for output for humans
				std::string type = "wrap_around";

				filename += type;
				filename += dottxt;

				vecut::write_into_file(filename, data); 
			}
			else {
				// output data in human readable format, along with frequency data

				// output the positive frequency components of the data set
				// set up the frequency space for the positive half of the FFT

				int N_fr = static_cast<int>(data.size()) / 4; // this may have been re-sized, so use this value instead of the one that was input
											// divide by 4 because array is of length 2*N
				std::vector<double> fr_vals(N_fr, 0.0);

				double delta_fr, fr0 = 0.0, fr_final = 1.0 / (2.0*pos_spac);

				delta_fr = (fr_final - fr0) / (static_cast<double>(N_fr - 1));

				/*std::cout << "\nFrequency Space\n";
				std::cout << "N = " << N_fr << ", delta-T = " << pos_spac << ", 1/2T = " << fr_final << "\n";
				std::cout << "f0 = " << fr0 << " , ff = " << fr_final << " , df = " << delta_fr << "\n";*/

				for (int i = 0; i < N_fr; i++) {
					fr_vals[i] = fr0;
					fr0 += delta_fr;
				}				

				std::string fr_file = filename + "_Frq_data" + dottxt; // filename for frequency data

				vecut::write_into_file(fr_file, fr_vals);

				// output the positive frequency components of the computed FFT

				std::string fft_file = filename + "_Abs_FFT_data" + dottxt; // filename for absolute value of FFT data

				std::ofstream write(fft_file, std::ios_base::out, std::ios_base::trunc);

				if (write.is_open()) {

					for (size_t i = 0; i < data.size() / 2; i += 2) {
						write << std::setprecision(10) << template_funcs::Pythag(data[i], data[i+1]) << "\n";
					}

					write.close();
				}

				fr_vals.clear();
			}
		}
		else {
			std::string reason; 
			reason = "Error: void fft::output_data(std::vector<double> data, double pos_spac, std::string &filename, bool wrap_around)\n";
			if (data.empty()) reason += "data is empty\n"; 
			if (filename == empty_str || !useful_funcs::valid_filename_length(filename)) reason += "Filename: " + filename + " is not valid\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

// Private member method definitions

void fft::format_four1(std::vector<double> &data)
{
	// format a real data set of length N for input into four1
	// this converts the real data set of length N into a complex array of length 2*N
	// store values in the form (real, imag, real, imag, real, imag, ....), with all imag = 0
	// no padding is performed

	try {
		if (data.size() > 1) {
			std::vector<double> tmp_data(2 * data.size(), 0.0);

			int j = 0;
			for (size_t i = 0; i < tmp_data.size(); i += 2) {
				tmp_data[i] = data[j];
				j++;
			}

			// Store tmp_data in data
			data = tmp_data;

			tmp_data.clear();
		}
		else {
			std::string reason;
			reason = "Error: void fft::format_four1(std::vector<double> &data)\n";
			reason += "Vector input size not valid\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fft::pad_data(std::vector<double> &data, unsigned long &nn, bool CMPLX_ARR)
{
	// Pad a data set with zeroes if the length of the data set is not a power of two
	// Only need to add zeroes to the data set if nn is not a power of two

	try {
		if (data.size() > 1) {
			if (!useful_funcs::is_POT(nn)) {

				// Convert nn to the next highest power of two
				nn = useful_funcs::next_POT(nn);

				std::vector<double> tmp_data;

				if (CMPLX_ARR) {
					// data is an array of complex numbers and must have length 2*nn
					tmp_data.resize(2 * nn, 0.0);
					//tmp_data.zero();

					for (size_t i = 0; i < data.size(); i++) {
						tmp_data[i] = data[i];
					}
				}
				else {
					// data is an array of real numbers and must have length nn
					tmp_data.resize(nn, 0.0);
					//tmp_data.zero();

					for (size_t i = 0; i < data.size(); i++) {
						tmp_data[i] = data[i];
					}
				}

				// Store tmp_data in data
				data = tmp_data;

				tmp_data.clear();
			}
		}
		else {
			std::string reason;
			reason = "Error: void fft::pad_data(std::vector<double> &data, unsigned long &nn, bool CMPLX_ARR)\n";
			reason += "Vector input size not valid\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fft::four1(std::vector<double> &data, unsigned long nn, int isign)
{
	// Replaces data[0..2*nn-1] by its discrete Fourier transform, if isign is input as 1
	// Replace data[0..2*nn-1] by nn times its inverse dicrete Fourier transform, if isign is input as -1
	// data is a complex array of length nn, or equivalently, a real array of length 2*nn
	// nn MUST be an integer power of 2

	unsigned long n, mmax, m, j, istep, i;
	// Need double precision for trigonometric recurrences
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	n = nn << 1;
	j = 1;
	for (i = 1; i < n; i += 2) {
		if (j > i) {
			// This is the bit reversal section of the routine
			template_funcs::SWAP(data[j-1], data[i-1]);
			template_funcs::SWAP(data[j], data[i]); // Exchange the two complex numbers
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while (n > mmax) { // Outer loop executed \log_{2}{nn} times
		istep = mmax << 1;
		//theta=isign*(6.28318530717959/mmax); // Initialise the trigonometric recurrences
		theta = isign * (Two_PI / mmax); // Initialise the trigonometric recurrences
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) { // Here are the two nested inner loops
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr * data[j-1] - wi * data[j]; // Danielson-Lanczos formula
				tempi = wr * data[j] + wi * data[j-1];
				data[j-1] = data[i-1] - tempr;
				data[j] = data[i] - tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wr = (wtemp = wr)*wpr - wi * wpi + wr; // trigonometric recurrence
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}
}