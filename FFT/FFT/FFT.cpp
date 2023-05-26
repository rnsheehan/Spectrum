#include "pch.h" // use stdafx.h in Visual Studio 2017 and earlier

double convert_dBm_to_mW(double dBm_val)
{
	// convert a dBm power reading to a mW power reading
	return pow(10.0, (dBm_val / 10.0));
}

double convert_mW_to_dBm(double mW_val)
{
	// convert a mW power reading to dBm power reading
	if (mW_val > 1.0e-9) {
		return 10.0 * log10(mW_val);
	}
	else {
		return -90.0;
	}
}

void compute_FFT(unsigned long& N_spctr_data, double& spctr_hor_spacing, std::vector<double>& spctr_data, int& N_fft_data, std::vector<double>& fft_data, std::vector<double>& fft_abcissae)
{
	// Compute the FFT of spctr_data
	// N_spctr_data is the number of measured spectral data points
	// spctr_hor_spacing is the spacing between the measured spectral data points on the horizontal axis
	// spctr_data is the measured spectral data set
	// N_fft_data will contain the number of computed FFT data points
	// fft_data will contain the positive frequency components of the computed FFT, phase information is not retained here
	// fft_abcissae will contain the horizontal positions for the fft_data in the transformed space
	// R. Sheehan 19 - 11 - 2021

	try {
		bool c1 = N_spctr_data > 0 ? true : false;
		bool c2 = spctr_hor_spacing > 0 ? true : false;
		bool c3 = spctr_data.size() == N_spctr_data ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {

			fft calc;

			calc.compute_transform(N_spctr_data, spctr_hor_spacing, spctr_data, N_fft_data, fft_data, fft_abcissae);
		}
		else {
			std::string reason = "Error: void testing::compute_FFT()\n";
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

void Spectrum_FFT(unsigned long& N_spctr_data, double spctr_hor_spacing, double spctr_data[], int& N_fft_data, double fft_data[], double fft_abcissae[])
{
	// Interface for accessing the FFT computation methods
	// N_spctr_data is the number of measured spectral data points
	// spctr_hor_spacing is the spacing between the measured spectral data points on the horizontal axis
	// spctr_data is the measured spectral data set
	// N_fft_data will contain the number of computed FFT data points
	// fft_data will contain all frequency components of the computed FFT
	// fft_abcissae will contain the horizontal positions for the fft_data in the transformed space
	// R. Sheehan 19 - 11 - 2021

	// Create std::vector for computing the fits
	std::vector<double> y(N_spctr_data, 0.0);

	// store the data in the vector containers instead of standard arrays
	// use vectors since that is what the FFT code uses
	double scale_fac = 1.0e+6;
	for (int i = 0; i < static_cast<int>(N_spctr_data); i++) {
		y[i] = scale_fac * convert_dBm_to_mW(spctr_data[i]); // convert from dBm scale to mW scale
	}

	// create variables for storing the computed results
	// use vectors since that is what the FFT code uses
	N_fft_data = 0;
	std::vector<double> fft_cont;
	std::vector<double> fft_x;

	// perform the FFT calculation
	compute_FFT(N_spctr_data, spctr_hor_spacing, y, N_fft_data, fft_cont, fft_x);

	// store the data in containers for output
	// use arrays since that is what LabVIEW uses
	// size of fft_data and fft_abcissae must be known a-priori
	// N_fft_data = 2 * NextPOT(N_spctr_data)

	// Arrays must have the correct size at calling time inside LabVIEW
	/*fft_data = new double[N_fft_data]; 
	fft_abcissae = new double[N_fft_data];*/

	for (int i = 0; i < N_fft_data; i++) {		
		fft_abcissae[i] = fft_x[i];
	}

	for (int i = 0; i < 2*N_fft_data; i++) {
		fft_data[i] = fft_cont[i]; 
	}
}