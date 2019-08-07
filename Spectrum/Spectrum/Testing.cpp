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

	calc.output_data(spctdata, delta_t, spctfile); 
}

void testing::laser_cavity_length()
{
	// Can you deduce the length of a laser cavity from the FFT of its spectrum
	// R. Sheehan 6 - 8 - 2019

	std::string timefile = "Actual_Frequency.txt";
	std::string spctfile = "Spectrum_I_150.txt";

	std::vector<double> timedata; int ntimes = 0;
	std::vector<double> spctdata; int nspct = 0;

	vecut::read_into_vector(timefile, timedata, ntimes);

	vecut::read_into_vector(spctfile, spctdata, nspct);

	double delta_t = (timedata[ntimes - 1] - timedata[0]) / (ntimes-1);

	unsigned long nn = nspct;

	fft calc;

	calc._four1(spctdata, nn); // compute the FFT of spctdata

	calc.output_data(spctdata, fabs(delta_t), spctfile);	
}