#ifndef CONVOLUTIONS_H
#define CONVOLUTIONS_H

// Implementation of an object used to compute the convolution and deconvolution of data sets
// Based on the functions given in NRinC, Ch 13
// R. Sheehan 23 - 3 - 2015

// convol_deconvol inherits fft which means that all public available to fft are available to convol_deconvol

class convol_deconvol : fft {
public:
	// Constructors
	convol_deconvol();

	//Methods
	// Add methods to check for errors etc. 

	void _convlv(std::vector<double>& data, unsigned long n, std::vector<double>& respns, unsigned long m, int isign, std::vector<double>& ans);

private:

	void convlv(std::vector<double>& data, unsigned long n, std::vector<double>& respns, unsigned long m, int isign, std::vector<double>& ans);
};

#endif
