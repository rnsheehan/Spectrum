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

	//std::string timefile = "Ed_Laser_WL.csv";
	//std::string spctfile = "Ed_Laser_Pow.csv";

	/*std::string timefile = "PM_Laser_WL.csv";
	std::string spctfile = "PM_Laser_Pow.csv";*/

	//std::string timefile = "WL_Meas_Shallow_I_90.txt";
	//std::string spctfile = "Pow_Meas_Shallow_I_90.txt";

	//std::string timefile = "Wavelength.txt";
	//std::string spctfile = "Spectrum_I_150.txt";

	std::string timefile = "20C_wave.dat";
	std::string spctfile = "50mA.dat";

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

	calc.output_data(spctdata, delta_t, spctfile, ".dat");	
}