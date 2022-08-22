#ifndef TESTING_H
#define TESTING_H

// Declaration of namespace that implements methods to test operation of FFT algorithm

namespace testing {
	void test_rotate(); 

	void sample_FFT_calculation(); 

	void laser_FFT(); 

	void laser_LLM(); 

	void compute_FFT(unsigned long& N_spctr_data, double& spctr_hor_spacing, std::vector<double>& spctr_data, int& N_fft_data, std::vector<double>& fft_data, std::vector<double>& fft_abcissae);

	void compute_FFT_test(); 

	void example_calculations();

	void inverse_FFT_test();

	void real_ft_test(); 

	void two_ft_test(); 

	void sine_wave(int Nsmpls, double Lt, double f1, double f2);

	void cosine_wave(int Nsmpls, double Lt, double f1, double f2); 

	void gaussian(int Nsmpls, double Lt, double mean, double stdev); 

	void exponential(int Nsmpls, double Lt, double centre, double decay);

	void tophat(int Nsmpls, double Lt, double centre, double width); 

	void tophat_alt(int Nsmpls, double Lt, double centre, double width);

	void heaviside(int Nsmpls, double Lt, double centre);
	
	void signum(int Nsmpls, double Lt, double centre);
}

#endif
