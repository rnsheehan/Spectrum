#ifndef ATTACH_H
#include "Attach.h"
#endif

// Constructors
convol_deconvol::convol_deconvol()
{
	// Default Constructor
}

void convol_deconvol::_convlv(std::vector<double>& data, unsigned long n, std::vector<double>& respns, unsigned long m, int isign, std::vector<double>& ans)
{
	// Calls convlv with exception handling

	// Convolves or deconvolves a real data set data[0..n-1] (including any user supplied padding)
	// with a response function respns[0..m-1]. 

	// The response function must be converted wrap-around order. This is done here before respns is passed to convlv
	// Wrap-around order means that the first half of the array contains the impulse reponse function at negative times, counting down 
	// from the highest element respns[m-1]. 

	// On input isign is +1 for convolution, -1 for deconvolution. 

	// The answer is returned in the first n components of ans[0..n-1]. 
	// n MUST be an integer power of 2.

	// "Whomsoever diggeth a pit, shall fall in it1
	// Whomsoever diggeth a pit, shall bury in it!
	// If you have a big tree, we have a small axe!"
	// "Small Axe", The Wailers
	// https://www.youtube.com/watch?v=WJJeLFvYsT0

	try {

		if (!useful_funcs::is_POT(n)) {
			bool CMPLX_ARR = false;
			unsigned long n_orig = n; // store the original length in case you need to resize each array
			unsigned long n2 = 2 * n; // store the value of the size of the array ans[0..2n-1]

			pad_data(data, n, CMPLX_ARR); // data[0..n-1]

			n = n_orig; // re-set the value of n in case it was altered by the last call to pad_data

			//pad_data(respns, n, CMPLX_ARR); // respns[0..n-1], the value stored in n may be altered after this call

			pad_data(ans, n, CMPLX_ARR); // ans[0..n-1], call to pad_data will change n to next POT
		}

		// Convert respns to wrap-around order here
		// Take M/2 values to the left for t < 0 and M/2 values to the right of mid_indx and store them in wrap around format
		//int count = 0;
		//int mid_indx = n / 2;
		//int m2 = m / 2;
		//for (int i = (m + 3) / 2; i <= m; i++) {
		//	respns[i] = respns[mid_indx - m2 + count];
		//	count++;
		//}

		//for (int i = 1; i < (m + 3) / 2; i++) {
		//	respns[i] = respns[mid_indx + i + 1]; // take data for t > 0
		//}

		// Convert respns to wrap-around order here
		vecut::wrap_around_conversion(respns); 

		// check the length of each array
		bool c1 = fabs(data.size() - static_cast<int>(n)) == 0;
		bool c2 = fabs(respns.size() - static_cast<int>(m)) == 0;
		bool c3 = fabs(ans.size() - static_cast<int>(n)) == 0;
		bool c4 = useful_funcs::is_POT(n); 
		bool c5 = fabs(isign) == 1; 
		bool c6 = m <= n; 
		bool c10 = c1 && c2 && c3 && c4 && c5 && c6; 

		if (c10) {

			convlv(data, n, respns, m, isign, ans);
		}
		else {
			std::string reason;
			reason = "Error: void convol_deconvol::_convlv(std::vector<double>& data, unsigned long n, std::vector<double>& respns, unsigned long m, int isign, std::vector<double>& ans)\n";
			if (!c1 || !c2 || !c3) {
				reason += "Input arrays do not have the correct dimensions\n";
				reason += "m = " + template_funcs::toString(m) + "\n";
				reason += "n = " + template_funcs::toString(n) + "\n";
				reason += "data.size() = " + template_funcs::toString(data.size()) + "\n";
				reason += "respns.size() = " + template_funcs::toString(respns.size()) + "\n";
				reason += "ans.size() = " + template_funcs::toString(ans.size()) + "\n";

			}
			if (!c4) {
				reason += "n is not a power of 2\n";
			}
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}

}

void convol_deconvol::convlv(std::vector<double>& data, unsigned long n, std::vector<double>& respns, unsigned long m, int isign, std::vector<double>& ans)
{
	// Convolves or deconvolves a real data set data[0..n-1] (including any user supplied padding)
	// with a response function respns[0..m-1]. 

	// The response function must be stored in wrap-around order in the first m elements of respns, where m is an odd integer <=n.
	// Wrap-around order means that the first half of the array contains the impulse reponse
	// function at negative times, counting down from the highest element respns[m]. 

	// On input isign is +1 for convolution, -1 for deconvolution. 

	// The answer is returned in the first n components of ans. 
	// However, ans must be supplied in the calling program with dimensions ans[1..2*n], 
	// for consistency with twofft	

	// n MUST be an integer power of 2. 

	// For response data in input
	/*cout<<"response data for times t > 0 contained in elements 1 <= i < "<<(m+3)/2<<"\n";
	cout<<"response data for times t <= 0 contained in elements "<<(m+3)/2<<" <= i <= "<<m<<"\n"; */

	// "I'm a well fed man so I grow up strong. 
	// I hope you overstand!"
	// "Favourite Dish", Lee 'Scratch' Perry
	// https://www.youtube.com/watch?v=mqH4Hoo011Y

	try {

		unsigned long i, no2;
		double dum, mag2;

		std::vector<double> temp(n, 0.0);
		temp[0] = respns[0];

		// put respns in array of length n
		for (i = 1; i < (m + 1) / 2; i++) {
			temp[i] = respns[i];
			temp[n - i] = respns[m - i];
		}

		// pad with zeros, don't need this since you've tmp is already filled
		/*for (i = (m + 1) / 2; i < n - (m - 1) / 2; i++) {
			temp[i] = 0.0;
		}*/

		for (i = 0; i < n; i++) {
			ans[i] = data[i];
		}

		// FFT both data sets
		_realft(ans, n, isign);
		_realft(temp, n, isign);

		no2 = n >> 1; // n02 = n/2

		if (isign == 1) { // multiply FFTs to convolve
			for (i = 2; i < n; i += 2) {
				dum = ans[i]; 
				ans[i] = (ans[i] * temp[i] - ans[i + 1] * temp[i + 1]) / no2; 
				ans[i + 1] = (ans[i+1] * temp[i] + dum * temp[i+1]) / no2; 
			}
			ans[0] = ans[0] * temp[0] / no2; 
			ans[1] = ans[1] * temp[1] / no2; 
		}
		else if (isign == -1) { // divide FFTs to deconvolve
			for (i = 2; i < n; i += 2) {
				if ((mag2 = template_funcs::DSQR(temp[i]) + template_funcs::DSQR(temp[i + 1])) == 0) {
					std::string reason;
					reason = "Error: void convol_deconvol::convlv(std::vector<double>& data, unsigned long n, std::vector<double>& respns, unsigned long m, int isign, std::vector<double>& ans)\n";
					reason += "Deconvolving at response zero in CONVLV\n";
					throw std::invalid_argument(reason);
				}
				else {
					dum = ans[i]; 
					ans[i] = (ans[i] * temp[i] + ans[i + 1] * temp[i + 1]) / mag2 / no2; 
					ans[i+1] = (ans[i+1] * temp[i] - dum * temp[i + 1]) / mag2 / no2; 
				}
			}
			if (temp[0] == 0 || temp[1] == 0.0) {
				std::string reason;
				reason = "Error: void convol_deconvol::convlv(std::vector<double>& data, unsigned long n, std::vector<double>& respns, unsigned long m, int isign, std::vector<double>& ans)\n";
				reason += "Deconvolving at response zero in CONVLV\n";
				throw std::invalid_argument(reason);
			}
			else {
				ans[0] = ans[0] / temp[0] / no2; 
				ans[1] = ans[1] / temp[1] / no2; 
			}
		}
		else {
			std::string reason;
			reason = "Error: void convol_deconvol::convlv(std::vector<double>& data, unsigned long n, std::vector<double>& respns, unsigned long m, int isign, std::vector<double>& ans)\n";
			reason += "No meaning for ISIGN in CONVLV\n";
			throw std::invalid_argument(reason);
		}

		// Inverse transform back to time domain
		_realft(ans, n, -1); 

		temp.clear();
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}