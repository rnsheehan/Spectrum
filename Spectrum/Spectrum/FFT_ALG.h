#ifndef FFT_ALG_H
#define FFT_ALG_H

// Implementation of an object that calls the FFT algorithms as given in NRinC, Ch. 12
// R. Sheehan 23 - 3 - 2015

class fft {
public:
	fft();

	void _four1(std::vector<double> &data, unsigned long &nn, int isign = 1, bool FORMAT_DATA = true);

	void _twofft(std::vector<double>& data1, std::vector<double>& data2, std::vector<double>& fft1, std::vector<double>& fft2, unsigned long n);

	void _realft(std::vector<double>& data, unsigned long n, int isign);

	void output_pos_data(std::vector<double> data, double smpl_spac, std::string &filename, std::string extension, bool output_type = STNDRD);

	void output_data(std::vector<double> data, double smpl_spac, std::string &filename, std::string extension, int isign = 1);

	void output_wrap_around(std::vector<double> data, std::string& filename, std::string extension); 

	void compute_transform(unsigned long& N_spctr_data, double& spctr_hor_spacing, std::vector<double>& spctr_data, int& N_fft_data, std::vector<double>& fft_data, std::vector<double>& fft_abcissae);

private:
	void four1(std::vector<double> &data, unsigned long nn, int isign);

	void twofft(std::vector<double> &data1, std::vector<double> &data2, std::vector<double> &fft1, std::vector<double> &fft2, unsigned long n);

	void realft(std::vector<double> &data, unsigned long n, int isign);

	void format_four1(std::vector<double> &data);

	void pad_data(std::vector<double> &data, unsigned long &nn, bool CMPLX_ARR);

	void create_freq_values(int & N_smpls, double &smpl_spac, std::vector<double> &fr_vals, int isign = 1, bool loud = false);

	void create_pos_freq_values(int &N_fr, double &smpl_spac, std::vector<double> &fr_vals, bool loud = false);

private:
	unsigned long arr_size;
};

#endif
