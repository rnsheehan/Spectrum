#ifndef FFT_ALG_H
#define FFT_ALG_H

// Implementation of an object that calls the FFT algorithms as given in NRinC, Ch. 12
// R. Sheehan 23 - 3 - 2015

class fft {
public:
	fft();

	void _four1(std::vector<double> &data, unsigned long &nn, int isign = 1, bool FORMAT_DATA = true);

	void output_pos_data(std::vector<double> data, double smpl_spac, std::string &filename, std::string extension, bool output_type = STNDRD);

	void output_data(std::vector<double> data, double smpl_spac, std::string &filename, std::string extension);

	void compute_transform(unsigned long& N_spctr_data, double& spctr_hor_spacing, std::vector<double>& spctr_data, int& N_fft_data, std::vector<double>& fft_data, std::vector<double>& fft_abcissae);

private:
	void four1(std::vector<double> &data, unsigned long nn, int isign);

	void format_four1(std::vector<double> &data);

	void pad_data(std::vector<double> &data, unsigned long &nn, bool CMPLX_ARR);

	void create_freq_values(int & N_smpls, double &smpl_spac, std::vector<double> &fr_vals, bool loud = false);

	void create_pos_freq_values(int &N_fr, double &smpl_spac, std::vector<double> &fr_vals, bool loud = false);

private:
	unsigned long arr_size;
};

#endif
