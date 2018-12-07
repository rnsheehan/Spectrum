#ifndef FFT_ALG_H
#define FFT_ALG_H

// Implementation of an object that calls the FFT algorithms as given in NRinC, Ch. 12
// R. Sheehan 23 - 3 - 2015

class fft {
public:
	fft();

	void _four1(std::vector<double> &data, unsigned long &nn, int isign, bool FORMAT_DATA = true);

private:
	void four1(std::vector<double> &data, unsigned long nn, int isign);

	void format_four1(std::vector<double> &data);

	void pad_data(std::vector<double> &data, unsigned long &nn, bool CMPLX_ARR);

private:
	unsigned long arr_size;
};

#endif
