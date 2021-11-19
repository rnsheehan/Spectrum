#ifndef TESTING_H
#define TESTING_H

// Declaration of namespace that implements methods to test operation of FFT algorithm

namespace testing {
	void sample_FFT_calculation(); 

	void laser_FFT(); 

	void laser_LLM(); 

	void compute_FFT(unsigned long& N_spctr_data, double& spctr_hor_spacing, std::vector<double>& spctr_data, int& N_fft_data, std::vector<double>& fft_data, std::vector<double>& fft_abcissae);

	void compute_FFT_test(); 
}

#endif
