#ifndef ATTACH_H
#include "Attach.h"
#endif // !ATTACH_H

// The aim of this project is to implement the code needed to compute the FFT of an arbitrary data set
// The data I have in mind is that obtained from a high speed oscilloscope
// Generally, the data will be of the form (time_value, voltage_value)
// The code for computing the FFT is taken from NRinC, this code adapts an existing code to run using std::vector<double>
// One of the tasks required is to determine the SNR of the FFT from the trace data
// R. Sheehan 5 - 12 - 2018

int main() 
{
	//testing::sample_FFT_calculation(); 

	//testing::laser_FFT();  

	//testing::laser_LLM(); 

	//testing::compute_FFT_test(); 

	//testing::example_calculations(); 

	testing::inverse_FFT_test(); 

	//testing::real_ft_test(); 

	//testing::two_ft_test(); 

	std::cout << "Press enter to close\n"; 
	std::cin.get(); 

	return 0;
}
