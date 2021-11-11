#ifndef ATTACH_H
#include "Attach.h"
#endif

void testing::sample_FFT_calculation()
{
	// sample FFT calculations to test operation of FFT object

	std::string timefile = "Time_Data.txt"; 
	std::string spctfile = "Spec_Data.txt"; 

	std::vector<double> timedata; int ntimes = 0; 
	std::vector<double> spctdata; int nspct = 0; 

	vecut::read_into_vector(timefile, timedata, ntimes);

	vecut::read_into_vector(spctfile, spctdata, nspct); 

	double delta_t = timedata[1] - timedata[0]; 

	unsigned long nn = nspct; 

	fft calc; 

	calc._four1(spctdata, nn); // compute the FFT of spctdata

	calc.output_data(spctdata, delta_t, spctfile, dottxt); 
}

void testing::laser_FFT()
{
	// Can you deduce the length of a laser cavity from the FFT of its spectrum
	// R. Sheehan 6 - 8 - 2019

	std::string dir_name = "c:\\users\\robert\\Research\\FFT_Data\\Examples_2019\\";

	useful_funcs::set_directory(dir_name);

	//useful_funcs::get_directory();

	/*std::string timefile = "Ed_Laser_WL.csv";
	std::string spctfile = "Ed_Laser_Pow.csv";*/

	/*std::string timefile = "PM_Laser_WL.csv";
	std::string spctfile = "PM_Laser_Pow.csv";*/

	std::string timefile = "HSSCP_Time_2.csv";
	std::string spctfile = "HSSCP_SPCT_2.csv";

	//std::string timefile = "WL_Meas_Shallow_I_90.txt";
	//std::string spctfile = "Pow_Meas_Shallow_I_90.txt";

	/*std::string timefile = "Wavelength.txt";
	std::string spctfile = "Spectrum_I_150.txt";*/

	//std::string timefile = "20C_wave.dat";
	//std::string spctfile = "50mA.dat";

	std::vector<double> timedata; int ntimes = 0;
	std::vector<double> spctdata; int nspct = 0;

	vecut::read_into_vector(timefile, timedata, ntimes);

	vecut::read_into_vector(spctfile, spctdata, nspct);

	double delta_t; 

	delta_t = (timedata[ntimes - 1] - timedata[0]) / (ntimes-1);

	//delta_t = timedata[1] - timedata[0]; 

	unsigned long nn = nspct;

	fft calc;

	calc._four1(spctdata, nn); // compute the FFT of spctdata

	calc.output_data(spctdata, delta_t, spctfile, dotcsv);	
}

void testing::laser_LLM()
{
	// compute the FFT of the measured laser linewidth spectral data
	// the FFT should give you the autocorrelation function
	// the rate of decay of the autocorrelation should give you the laser coherence time
	// R. Sheehan 11 - 11 - 2021

	// Read in the measured spectral data
	std::string filename = "Smpl_LLM_1.txt";
	std::string spctfile = "LLM_FFT.txt";

	int n_rows, n_cols; 
	std::vector<std::vector<double>> the_data;

	vecut::read_into_matrix(filename, the_data, n_rows, n_cols, true);

	// convert the spectral data from dBm to mW scale and rescale it
	double scale_fac = 1.0e+6;  double f_start = 60.0; double f_end = 100.0;
	std::vector<double> xdata; 
	std::vector<double> ydata; 
	for (int i = 0; i < n_rows; i++) {
		if (the_data[i][0] > f_start && the_data[i][0] < f_end) {
			xdata.push_back(the_data[i][0]);
			ydata.push_back(scale_fac * pow(10.0, the_data[i][1] / 10.0)); // convert the spectral data from dBm to mW scale and rescale it
		}
	}

	double delta_t = xdata[1] - xdata[0]; // frequency spacing in units of GHz

	unsigned long nn = ydata.size();

	fft calc;

	calc._four1(ydata, nn); // compute the FFT of spctdata

	calc.output_data(ydata, delta_t, spctfile, dottxt);
}