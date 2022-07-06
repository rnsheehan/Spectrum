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
	// data is the array whose FFT will be computed
	// nn is the length of the original data set, it gets changed by the pad data algorithm if necessary
	// isign is a switch that tells the code to compute the FFT or inverse FFT
	// FORMAT_DATA is a boolean switch that tells the code to decide if data must be formatted for FFT
	// see format_four1 for details 	

	// Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1
	// Replace data[1..2*nn] by nn times its inverse dicrete Fourier transform, if isign is input as -1
	// data is a complex array of length nn, or equivalently, a real array of length 2*nn
	// nn MUST be an integer power of 2

	// data must be in the form output by format_four1 before fft algorithm can be applied
	// padding of data can be applied after the fact

	try {
		
		if ( !data.empty() ) {

			if (FORMAT_DATA) format_four1(data);

			bool CMPLX_ARR = true; 

			pad_data(data, nn, CMPLX_ARR); // data input into four1 must be in cmplx format

			if (useful_funcs::is_POT(nn) && abs(isign) == 1) {

				// nn must be a power of two
				// data must be in the form (real, imag, real, imag, real, imag, ....), with all imag = 0, hence length 2*nn
				// isign = 1 for FT
				four1(data, nn, isign);

				// scale by 1/N in the case of an inverse FT
				if(isign == -1) for(size_t i = 0; i < data.size(); i++) data[i] /= nn;
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

void fft::output_data(std::vector<double> data, double smpl_spac, std::string& filename, std::string extension, int isign)
{
	// Output the full FFT spectrum at positive and negative frequency components
	// R. Sheehan 4 - 7 - 2022

	try {
	
		if (smpl_spac > 0.0 && !data.empty() && filename != empty_str && useful_funcs::valid_filename_length(filename))	{

			// remove .txt from end of filename
			useful_funcs::remove_substring(filename, extension);

			int N_smpls = data.size() / 2; // No. of complex valued data points in the FFT

			std::vector<double> fr_vals(N_smpls, 0.0);

			create_freq_values(N_smpls, smpl_spac, fr_vals, isign, true); // compute the frequency space values

			// filename for frequency data
			std::string fr_file = filename + "_Frq_data" + extension;

			vecut::write_into_file(fr_file, fr_vals);

			// output the positive + negative frequency components of the computed FFT

			std::string fft_file = filename + "_FFT_data" + extension; // filename for FFT data

			std::ofstream write(fft_file, std::ios_base::out, std::ios_base::trunc);

			if (write.is_open()) {
				
				if (isign == 1) {
					// output negative frequency FFT data in the form real, imag, abs-value, phase
					for (size_t i = N_smpls; i < data.size() - 1; i += 2) {
						write << std::setprecision(10) << data[i] << " , " << data[i + 1] << " , " << template_funcs::Pythag(data[i], data[i + 1]) << " , " << atan2(data[i + 1], data[i]) << "\n";
					}
				}

				// output positive frequency FFT data in the form real, imag, abs-value, phase
				for (size_t i = 0; i < N_smpls; i += 2) {
					write << std::setprecision(10) << data[i] << " , " << data[i + 1] << " , " << template_funcs::Pythag(data[i], data[i + 1]) << " , " << atan2(data[i + 1], data[i]) << "\n";
				}

				if (isign == -1) {
					// output the "negative" frequency components for the inverse FT
					for (size_t i = N_smpls; i < data.size() - 1; i += 2) {
						write << std::setprecision(10) << data[i] << " , " << data[i + 1] << " , " << template_funcs::Pythag(data[i], data[i + 1]) << " , " << atan2(data[i + 1], data[i]) << "\n";
					}
				}

				write.close();
			}

			fr_vals.clear();

		}
		else {
			std::string reason;
			reason = "Error: void fft::output_data(std::vector<double> data, double smpl_spac, std::string& filename, std::string extension)\n";
			if (data.empty()) reason += "data is empty\n";
			if (filename == empty_str || !useful_funcs::valid_filename_length(filename)) reason += "Filename: " + filename + " is not valid\n";
			if (!(smpl_spac > 0.0)) reason += "position spacing is not positive\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fft::output_wrap_around(std::vector<double> data, std::string& filename, std::string extension)
{
	// Output the full FFT spectrum at positive and negative frequency components in its native wrap-around format
	// R. Sheehan 5 - 7 - 2022

	try {

		if (!data.empty() && filename != empty_str && useful_funcs::valid_filename_length(filename)) {

			// remove .txt from end of filename
			useful_funcs::remove_substring(filename, extension);

			// output the positive + negative frequency components of the computed FFT

			std::string fft_file = filename + "_FFT_data_wrap_around" + extension; // filename for FFT data

			vecut::write_into_file(fft_file, data); 
		}
		else {
			std::string reason;
			reason = "Error: void fft::output_wrap_around(std::vector<double> data, std::string& filename, std::string extension)\n";
			if (data.empty()) reason += "data is empty\n";
			if (filename == empty_str || !useful_funcs::valid_filename_length(filename)) reason += "Filename: " + filename + " is not valid\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fft::output_pos_data(std::vector<double> data, double smpl_spac, std::string &filename, std::string extension, bool output_type)
{
	// output the computed FFT data to a file
	// How do you handle the position data if the size of the data has been expanded? 
	// Don't need it because you are working in frequency space. This means that you should output frequency data
	// for the human readable format. 

	// Need two options for output
	// Standard option is to output directly from the FFT algorithm, this outputs in wrap-around order
	// Also want an option to output data in a human readable / plottable format

	// this assumes that FFT of data has been computed

	// Explanation of wrap-around order from NRinC
	// Input array data contains N complex samples in a real array of length 2N (N must be POT) with real and imag parts alternating
	// On output data contains complex Fourier spectrum at N values of frequency with real and imag parts alternating 
	// N.B. Array starts at zero frequency, steps up to most positive frequency
	// N.B. negative frequencies follow from second most negative frequency up to frequency just below zero

	try {
		if (smpl_spac > 0.0 && !data.empty() && filename != empty_str && useful_funcs::valid_filename_length(filename)) {
		
			// remove .txt from end of filename
			useful_funcs::remove_substring(filename, extension);

			if (output_type) {
				// output data in standard wrap-around order
				// never really want to use this order for output for humans
				std::string type = "_FFT_wrap_around";

				filename += type;
				filename += extension;

				vecut::write_into_file(filename, data); 
			}
			else {
				// output data in human readable format, along with frequency data

				// output the positive frequency components of the data set
				// set up the frequency space for the positive half of the FFT
				// 
				// Explanation
				// Suppose data.size() = 4096 => 2048 complex-valued data points
				// which gives 1024 positive frequency component and 1024 negative frequency components
				// N_pos_fr = data.size() / 4 = 1024 positive ( negative )  frequency components
				// The frequency space spans -(1/(2smpl_spac)) < f < +(1/(2smpl_spac)) or [-f_end, 0) + [0, f_end], f_end = (1/(2smpl_spac))
				// Frequency spacing in positive frequency space must be delta_fr = f_end / N_pos_fr = 1/(2 N_pos_fr smpl_spac)
				// Frequency spacing over whole range will be given by delta_fr = 2 f_end / N_tot = 1/(N_tot smpl_spac) = 1/(2 N_pos_fr smpl_spac) since N_tot = 2 N_pos_fr

				int N_pos_fr = data.size() / 4; // this may have been re-sized, so use this value instead of the one that was input
				std::vector<double> fr_vals(N_pos_fr, 0.0);

				create_pos_freq_values(N_pos_fr, smpl_spac, fr_vals, true); // compute the frequency space values

				// filename for frequency data
				std::string fr_file = filename + "_Frq_data" + extension;

				vecut::write_into_file(fr_file, fr_vals);

				// output the positive frequency components of the computed FFT

				std::string fft_file = filename + "_FFT_data" + extension; // filename for FFT data

				std::ofstream write(fft_file, std::ios_base::out, std::ios_base::trunc);

				if (write.is_open()) {

					// output positive frequency FFT data in the form real, imag, abs-value, phase
					for (size_t i = 0; i < data.size() / 2; i += 2) {
						write << std::setprecision(10) << data[i] << " , "<< data[i+1]<<" , " << template_funcs::Pythag(data[i], data[i + 1]) << " , " << atan2(data[i + 1], data[i]) << "\n";
					}

					write.close();
				}

				fr_vals.clear();
			}
		}
		else {
			std::string reason; 
			reason = "Error: void fft::output_pos_data(std::vector<double> data, double smpl_spac, std::string &filename, bool wrap_around)\n";
			if (data.empty()) reason += "data is empty\n"; 
			if (filename == empty_str || !useful_funcs::valid_filename_length(filename)) reason += "Filename: " + filename + " is not valid\n"; 
			if (!(smpl_spac > 0.0)) reason += "position spacing is not positive\n"; 
			throw std::invalid_argument(reason); 
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fft::compute_transform(unsigned long& N_spctr_data, double& spctr_hor_spacing, std::vector<double>& spctr_data, int& N_fft_data, std::vector<double>& fft_data, std::vector<double>& fft_abcissae)
{
	// Compute the FFT of data
	// N_spctr_data is the number of measured spectral data points
	// spctr_hor_spacing is the spacing between the measured spectral data points on the horizontal axis
	// spctr_data is the measured spectral data set
	// Store the computed transform spectrum in {fft_abcissae, fft_data}
	// N_fft_data will contain the number of computed FFT data points
	// fft_data will contain the positive frequency components of the computed FFT, phase information is not retained here
	// fft_abcissae will contain the horizontal positions for the fft_data in the transformed space

	try {
		bool c1 = N_spctr_data > 0 ? true : false;
		bool c2 = spctr_hor_spacing > 0 ? true : false;
		bool c3 = spctr_data.size() == N_spctr_data ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			
			_four1(spctr_data, N_spctr_data); // compute the FFT, data re-formatting performed inside _four1

			N_fft_data = static_cast<int>(spctr_data.size() / 4);	// this may have been re-sized, so use this value instead of the one that was input
																	// divide by 4 because array is now of length 2*N
			
			fft_abcissae.resize(N_fft_data, 0.0); //resize the array to hold the transform space horizontal values

			fft_data.resize(N_fft_data, 0.0); //resize the array to hold the transform space horizontal values

			create_freq_values(N_fft_data, spctr_hor_spacing, fft_abcissae); // compute the transform space horizontal values

			// Store the absolute value of FFT spectrum in fft_data
			int count = 0; 
			for (size_t i = 0; i < spctr_data.size() / 2; i += 2) {
				fft_data[count] = template_funcs::Pythag(spctr_data[i], spctr_data[i + 1]);
				count++; 
			}
		}
		else {
			std::string reason = "Error: void fft::compute_transform()\n";
			if (!c1) reason += "N_spctr_data input is not correct\n";
			if (!c2) reason += "spctr_hor_spacing input is not correct\n";
			if (!c3) reason += "spctr_data is not correct size\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
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

			// test the array size to see if it is a POT
			if (!useful_funcs::is_POT(nn)) {
				// data has length nn which is not a POT this must be adjusted

				// Convert nn to the next highest power of two
				nn = useful_funcs::next_POT(nn);

				std::vector<double> tmp_data;

				if (CMPLX_ARR) {
					// data, which has already been converted to CMPLX_ARR, must be stored in an array of size 2nn, nn is a POT
					tmp_data.resize(2 * nn, 0.0);
					
					for (size_t i = 0; i < data.size(); i++) {
						tmp_data[i] = data[i];
					}
				}
				else {
					// data is an array of real numbers and must have length nn

					// what's the point of this option? is this to be used on the case of realft? R. Sheehan 5 - 7 - 2022

					tmp_data.resize(nn, 0.0);
					
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

	// Explanation of wrap-around order from NRinC
	// Input array data contains N complex samples in a real array of length 2N (N must be POT) with real and imag parts alternating
	// On output data contains complex Fourier spectrum at N values of frequency with real and imag parts alternating 
	// N.B. Array starts at zero frequency, steps up to most positive frequency
	// N. B. negative frequencies follow from second most negative frequency up to frequency just below zero

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
		m = n >> 1; // n>>1 => n -> n/2
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1; // n>>1 => n -> n/2
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

void fft::twofft(std::vector<double>& data1, std::vector<double>& data2, std::vector<double>& fft1, std::vector<double>& fft2, unsigned long n)
{
	// Given two real input arrays data1[1..n] and data2[1..n], this routine calls four1
	// ans returns two complex output arrays, fft1[1..2n] and fft2[1..2n], each of complex
	// length n, i.e. real length 2*n, which contain the discrete Fourier Transforms of the respective data
	// n MUST be an integer power of 2

	unsigned long nn3, nn2, jj, j;
	double rep, rem, aip, aim;

	nn3 = 1 + (nn2 = 2 + n + n);
	for (j = 0, jj = 1; j < n; j++, jj += 2) {
		fft1[jj - 1] = data1[j];
		fft1[jj] = data2[j];
	}
	four1(fft1, n, 1);
	fft2[0] = fft1[1];
	fft1[1] = fft2[1] = 0.0;
	for (j = 2; j < n + 1; j += 2) {
		rep = 0.5 * (fft1[j] + fft1[nn2 - j]);
		rem = 0.5 * (fft1[j] - fft1[nn2 - j]);
		aip = 0.5 * (fft1[j + 1] + fft1[nn3 - j]);
		aim = 0.5 * (fft1[j + 1] - fft1[nn3 - j]);
		fft1[j] = rep;
		fft1[j + 1] = aim;
		fft1[nn2 - j] = rep;
		fft1[nn3 - j] = -aim;
		fft2[j] = aip;
		fft2[j + 1] = -rem;
		fft2[nn2 - j] = aip;
		fft2[nn3 - j] = rem;
	}
}

void fft::realft(std::vector<double>& data, unsigned long n, int isign)
{
	// Calculates the Fourier transform of a set of n real-valued data points
	// Replaces this data, stored in data[1..n], by the positive frequency half of its complex Fourier transform
	// The real valued first and last components of the complex transform are returned as elements
	// data[1] and data[2] respectively.
	// This routine also calculates the inverse transform of a complex data array if it is the transform of real 
	// data, result in this case must be multiplied by 2/n
	// n MUST be an integer power of 2

	unsigned long i, i1, i2, i3, i4, np3;
	double c1 = 0.5, c2, h1r, h1i, h2r, h2i;
	double wr, wi, wpr, wpi, wtemp, theta;

	//theta=3.141592653589793/(double) (n>>1);
	theta = PI / (double)(n >> 1); // n>>1 => n -> n/2
	if (isign == 1) {
		c2 = -0.5;
		four1(data, n >> 1, 1); // n>>1 => n -> n/2
	}
	else {
		c2 = 0.5;
		theta = -theta;
	}
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0 + wpr;
	wi = wpi;
	np3 = n + 3;
	for (i = 1; i < (n >> 2); i++) { // n>>2 => n -> n/4
		i4 = 1 + (i3 = np3 - (i2 = 1 + (i1 = i + i - 1)));
		h1r = c1 * (data[i1] + data[i3]);
		h1i = c1 * (data[i2] - data[i4]);
		h2r = -c2 * (data[i2] + data[i4]);
		h2i = c2 * (data[i1] - data[i3]);
		data[i1] = h1r + wr * h2r - wi * h2i;
		data[i2] = h1i + wr * h2i + wi * h2r;
		data[i3] = h1r - wr * h2r + wi * h2i;
		data[i4] = -h1i + wr * h2i + wi * h2r;
		wr = (wtemp = wr) * wpr - wi * wpi + wr;
		wi = wi * wpr + wtemp * wpi + wi;
	}
	if (isign == 1) {
		data[0] = (h1r = data[0]) + data[1];
		data[1] = h1r - data[1];
	}
	else {
		data[0] = c1 * ((h1r = data[0]) + data[1]);
		data[1] = c1 * (h1r - data[1]);
		four1(data, n >> 1, -1); // n>>1 => n -> n/2
	}
}

void fft::create_pos_freq_values(int &N_pos_fr, double &smpl_spac, std::vector<double> &fr_vals, bool loud)
{
	// Create the set of positive frequency values that correspond to a straightforward FFT calculation
	// R. Sheehan 8 - 8 - 2019

	// output the positive frequency components of the data set
	// set up the frequency space for the positive half of the FFT
	
	// Explanation
	// Suppose data.size() = 4096 => 2048 complex-valued data points
	// which gives 1024 positive frequency component and 1024 negative frequency components
	// N_pos_fr = data.size() / 4 = 1024 positive ( negative )  frequency components
	// The frequency space spans -(1/(2smpl_spac)) < f < +(1/(2smpl_spac)) or [-f_end, 0) + [0, f_end], f_end = (1/(2 smpl_spac))
	// Frequency spacing in positive frequency space must be delta_fr = f_end / N_pos_fr = 1/(2 N_pos_fr smpl_spac)
	// Frequency spacing over whole range will be given by delta_fr = 2 f_end / N_tot = 1/(N_tot smpl_spac) = 1/(2 N_pos_fr smpl_spac) since N_tot = 2 N_pos_fr

	try {
		if (smpl_spac > 0.0 && static_cast<int>(fr_vals.size()) == N_pos_fr) {
			
			// I don't think this is correct
			//double delta_fr, fr0 = 0.0, fr_final = 1.0 / (2.0*smpl_spac);
			//delta_fr = (fr_final - fr0) / (static_cast<double>(N_fr - 1));

			// this is the definition supplied in NRinC
			double fr0 = 0.0, delta_fr = 1.0 / (2.0 * N_pos_fr * smpl_spac), fr_final = 1.0/(2.0*smpl_spac);
			//delta_fr = (fr_final - fr0) / (static_cast<double>(N_pos_fr));

			if(loud){
				std::cout << "\nPositive Frequency Space\n";
				std::cout << "No. positive freq. cpts N_{pos} = " << N_pos_fr << ", delta-T = " << smpl_spac << ", 1/2T = " << fr_final << "\n";
				std::cout << "f0 = " << fr0 << " , ff = " << fr_final << " , df = " << delta_fr << "\n\n"; 
			}			

			for (int i = 0; i < N_pos_fr; i++) {
				fr_vals[i] = fr0;
				fr0 += delta_fr;
			}
		}
		else {
			std::string reason;
			reason = "Error: void fft::create_pos_freq_values(int &N_fr, double &smpl_spac, std::vector<double> &fr_vals, bool loud)\n";
			reason += "smpl_spac or vector input size not valid\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument &e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void fft::create_freq_values(int& N_smpls, double& smpl_spac, std::vector<double>& fr_vals, int isign, bool loud)
{
	// create the full list of positive and negative frequency values in the case of FFT
	// create the full list of time values in the case of IFT
	// R. Sheehan 4 - 7 - 2022

	try {
	
		if (smpl_spac > 0.0 && static_cast<int>(fr_vals.size()) == N_smpls && abs(isign) == 1)	{
			// Compute the sample spacing in frequency space
			// delta_f = 1 / (N_smpls * delta), delta = sample spacing in time space
			double delta_fr = 1.0 / (N_smpls * smpl_spac), fr_final = 0.0, fr0 = 0.0;

			if (isign == 1) {
				fr_final = 1.0 / (2.0 * smpl_spac);  fr0 = -fr_final;
			}

			if (isign == -1) {
				fr_final = N_smpls * delta_fr;  fr0 = 0.0;
			}

			if (loud && isign == 1) {
				std::cout << "\nFrequency Space\n";
				std::cout << "No. freq. cpts N = " << N_smpls << ", delta-T = " << smpl_spac << ", 1/2T = " << fr_final << "\n";
				std::cout << "f0 = " << fr0 << " , ff = " << fr_final << " , df = " << delta_fr << "\n\n";
			}

			if (loud && isign == -1) {
				std::cout << "\nTime Space\n";
				std::cout << "No. time cpts N = " << N_smpls << ", delta-f = " << smpl_spac << "\n";
				std::cout << "t0 = " << fr0 << " , tf = " << fr_final << " , df = " << delta_fr << "\n\n";
			}

			for (int i = 0; i < N_smpls; i++) {
				fr_vals[i] = fr0;
				fr0 += delta_fr;
			}
		}
		else {
			std::string reason;
			reason = "Error: void fft::create_freq_values(int &N_smpls, double &smpl_spac, std::vector<double> &fr_vals, bool loud)\n";
			reason += "smpl_spac or vector input size not valid\n";
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}