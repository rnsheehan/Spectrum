#ifndef ATTACH_H
#include "Attach.h"
#endif

void testing::test_rotate()
{
	// convert the elements in a vector from wrap-around order to standard order
	// take elements from positions [size()/2, size()-1] and move them to positions [0, -1 + size()/2]
	// take elements from positions [0, -1 + size()/2] and move them to positions [size()/2, size()-1]
	// It would be nice if you could do this in place and it turns out that you can
	// While looking for a generic function to convert the output of the FFT algorithm from wrap-around to standard ordering
	// I encountered the std::rotate method
	// see https://cplusplus.com/reference/algorithm/rotate/ for details
	// R. Sheehan 22 - 8 - 2022

	bool USE_CUSTOM_ROTATE = true; 

	std::vector<double> myvector;

	// set some values:
	for (int i = 1; i < 14; ++i) myvector.push_back(i); // 1 2 3 4 5 6 7 8 9 "standard order"

	std::cout << "myvector initially contains:";
	for (std::vector<double>::iterator it = myvector.begin(); it != myvector.end(); ++it)std::cout << ' ' << *it;
	std::cout << '\n';

	if (USE_CUSTOM_ROTATE) {
		vecut::wrap_around_conversion(myvector, false); 
	}
	else {
		std::rotate(myvector.begin(), myvector.begin() + myvector.size() / 2, myvector.end()); // 5 6 7 8 9 1 2 3 4 "wrap-around order"
	}

	// print out content:
	std::cout << "myvector after rotation contains:";
	for (std::vector<double>::iterator it = myvector.begin(); it != myvector.end(); ++it)std::cout << ' ' << *it;
	std::cout << '\n';

	if (USE_CUSTOM_ROTATE) {
		vecut::wrap_around_conversion(myvector, true);
	}
	else {
		if ((myvector.size()) % 2 == 0) {
			std::cout << (myvector.size()) % 2 << "\n";
			std::rotate(myvector.begin(), myvector.begin() + myvector.size() / 2, myvector.end()); // 1 2 3 4 5 6 7 8 9 "standard order"
		}
		else {
			std::rotate(myvector.begin(), myvector.begin() + (myvector.size() / 2) + 1, myvector.end()); // 1 2 3 4 5 6 7 8 9 "standard order"
		}
	}

	std::cout << "myvector after rotation contains:";
	for (std::vector<double>::iterator it = myvector.begin(); it != myvector.end(); ++it)std::cout << ' ' << *it;
	std::cout << '\n';
}

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

	/*std::string timefile = "Ed_Laser_WL.csv";
	std::string spctfile = "Ed_Laser_Pow.csv";*/

	/*std::string timefile = "PM_Laser_WL.csv";
	std::string spctfile = "PM_Laser_Pow.csv";*/

	std::string timefile = "HSSCP_Time_2.csv";
	std::string spctfile = "HSSCP_SPCT_2.csv";

	//std::string timefile = "WL_Meas_Shallow_I_90.txt";
	//std::string spctfile = "Pow_Meas_Shallow_I_90.txt";

	/*std::string timefile = "Wavelength.txt";
	std::string spctfile = "Spectrum_I_150.txt";*/

	//std::string timefile = "20C_wave.dat";
	//std::string spctfile = "50mA.dat";

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

	calc.output_data(spctdata, delta_t, spctfile, dotcsv);	
}

void testing::laser_LLM()
{
	// compute the FFT of the measured laser linewidth spectral data
	// the FFT should give you the autocorrelation function
	// the rate of decay of the autocorrelation should give you the laser coherence time
	// R. Sheehan 11 - 11 - 2021

	// Read in the measured spectral data
	std::string filename = "Smpl_LLM_1.txt";
	std::string spctfile = "LLM_FFT.txt";

	int n_rows, n_cols; 
	std::vector<std::vector<double>> the_data;

	vecut::read_into_matrix(filename, the_data, n_rows, n_cols, true);

	// convert the spectral data from dBm to mW scale and rescale it
	double scale_fac = 1.0e+6;  double f_start = 60.0; double f_end = 100.0;
	std::vector<double> xdata; 
	std::vector<double> ydata; 
	for (int i = 0; i < n_rows; i++) {
		if (the_data[i][0] > f_start && the_data[i][0] < f_end) {
			xdata.push_back(the_data[i][0]);
			ydata.push_back(scale_fac * pow(10.0, the_data[i][1] / 10.0)); // convert the spectral data from dBm to mW scale and rescale it
		}
	}

	double delta_t = xdata[1] - xdata[0]; // frequency spacing in units of MHz

	unsigned long nn = ydata.size();

	fft calc;

	calc._four1(ydata, nn); // compute the FFT of spctdata

	calc.output_data(ydata, delta_t, spctfile, dottxt);
}

void testing::compute_FFT(unsigned long& N_spctr_data, double& spctr_hor_spacing, std::vector<double>& spctr_data, int& N_fft_data, std::vector<double>& fft_data, std::vector<double>& fft_abcissae)
{
	// Compute the FFT of spctr_data
	// N_spctr_data is the number of measured spectral data points
	// spctr_hor_spacing is the spacing between the measured spectral data points on the horizontal axis
	// spctr_data is the measured spectral data set
	// N_fft_data will contain the number of computed FFT data points
	// fft_data will contain the positive frequency components of the computed FFT, phase information is not retained here
	// fft_abcissae will contain the horizontal positions for the fft_data in the transformed space
	// R. Sheehan 19 - 11 - 2021

	try {
		bool c1 = N_spctr_data > 0 ? true : false; 
		bool c2 = spctr_hor_spacing > 0 ? true : false; 
		bool c3 = spctr_data.size() == N_spctr_data ? true : false;
		bool c10 = c1 && c2 && c3; 

		if (c10) {
			
			fft calc;

			calc.compute_transform(N_spctr_data, spctr_hor_spacing, spctr_data, N_fft_data, fft_data, fft_abcissae); 
		}
		else {
			std::string reason = "Error: void testing::compute_FFT()\n"; 
			if (!c1) reason += "N_spctr_data input is not correct\n"; 
			if (!c2) reason += "spctr_hor_spacing input is not correct\n"; 
			if (!c3) reason += "spctr_data is not correct size\n"; 
			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		std::cerr << e.what();
	}
}

void testing::compute_FFT_test()
{
	// test the operation of compute_FFT
	// R. Sheehan 19 - 11 - 2021

	// Read in the measured spectral data
	std::string filename = "Smpl_LLM_1.txt";

	int n_rows, n_cols;
	std::vector<std::vector<double>> the_data;

	vecut::read_into_matrix(filename, the_data, n_rows, n_cols, true);

	// convert the spectral data from dBm to mW scale and rescale it
	double scale_fac = 1.0e+6;  double f_start = 60.0; double f_end = 100.0;
	std::vector<double> xdata;
	std::vector<double> ydata;
	for (int i = 0; i < n_rows; i++) {
		if (the_data[i][0] > f_start && the_data[i][0] < f_end) {
			xdata.push_back(the_data[i][0]);
			ydata.push_back(scale_fac * pow(10.0, the_data[i][1] / 10.0)); // convert the spectral data from dBm to mW scale and rescale it
		}
	}

	double spctr_hor_spacing = xdata[1] - xdata[0]; // frequency spacing in units of MHz

	unsigned long N_spctr_data = ydata.size();

	int N_fft_data = 0; 
	std::vector<double> fft_data; 
	std::vector<double> fft_abcissae;

	compute_FFT(N_spctr_data, spctr_hor_spacing, ydata, N_fft_data, fft_data, fft_abcissae); 

	std::string xfile = "fft_abscissae.txt";
	std::string yfile = "fft_ordinates.txt";

	vecut::write_into_file(xfile, fft_abcissae); 
	vecut::write_into_file(yfile, fft_data);
}

void testing::example_calculations()
{
	// perform some sample calculations
	// want to test that the numerical FFT computes the theoretical FFT
	// R. Sheehan 4 - 7 - 2022

	try {
		std::string the_dir = "c:\\users\\robertsheehan\\Research\\Notes\\FFT\\Examples";
		useful_funcs::set_directory(the_dir);

		// sine wave example
		int Nsmpls = 4096; 
		double Lx = 10, f1 = 12, f2 = 33; 
		
		//testing::sine_wave(Nsmpls, Lx, f1, f2);

		testing::gaussian(Nsmpls, Lx, 5, 0.5); 

		std::string func_str = "Gaussian"; // Unequal_Ampl_Sine_Wave 
		std::string timefile = func_str + "_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave time data
		std::string spctfile = func_str + "_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave voltage data

		std::vector<double> timedata; int ntimes = 0;
		std::vector<double> spctdata; int nspct = 0;

		vecut::read_into_vector(timefile, timedata, ntimes);

		vecut::read_into_vector(spctfile, spctdata, nspct);

		double delta_t = timedata[1] - timedata[0];

		unsigned long nn = nspct;

		fft calc;

		calc._four1(spctdata, nn); // compute the FFT of spctdata

		calc.output_data(spctdata, delta_t, spctfile, dottxt);

		spctdata.clear(); timedata.clear();

		// now perform the same calculation using realft and see what the output looks like
		//func_str = "Unequal_Ampl_Sine_Wave";
		//timefile = func_str + "_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave time data
		//spctfile = func_str + "_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave voltage data
		//vecut::read_into_vector(timefile, timedata, ntimes);
		//vecut::read_into_vector(spctfile, spctdata, nspct);
		//delta_t = timedata[1] - timedata[0];
		//nn = nspct;
		//calc._realft(spctdata, nn, true); 
		//calc.output_pos_data(spctdata, delta_t, spctfile, dottxt); 
		//spctdata.clear(); timedata.clear();
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void testing::inverse_FFT_test()
{
	// perform some sample calculations
	// want to test that the numerical FFT computes the theoretical FFT
	// R. Sheehan 4 - 7 - 2022

	try {
		std::string the_dir = "c:\\users\\robertsheehan\\Research\\Notes\\FFT\\Examples";
		useful_funcs::set_directory(the_dir);

		// sine wave example
		int Nsmpls = 4096;
		
		// I think the problems with the Signum and the Top-Hat functions are to do with the no. samples and the padding

		std::string func_str = "Gaussian";
		std::string timefile = func_str + "_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data
		std::string spctfile = func_str + "_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data

		std::vector<double> timedata; int ntimes = 0;
		std::vector<double> spctdata; int nspct = 0;

		vecut::read_into_vector(timefile, timedata, ntimes);

		vecut::read_into_vector(spctfile, spctdata, nspct);

		double delta_t = timedata[1] - timedata[0];

		unsigned long nn = nspct;

		fft calc;

		int isign = +1;
		bool FMT_DATA = true;

		calc._four1(spctdata, nn, isign, FMT_DATA); // compute the FFT of spctdata

		calc.output_data(spctdata, delta_t, spctfile, dottxt, isign); // output the data in appropriate format

		// Read the computed FFT data back into memory so that inverse FFT can be computed
		
		std::string freq_file = func_str + "_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + "_Frq_data" + dottxt;
		std::string fft_file = func_str + "_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + "_FFT_data" + dottxt;
		std::string fft_wrap = func_str + "_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + "_FFT_data_wrap_around" + dottxt;
		std::string ift_file = func_str + "_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + "_IFT" + dottxt;

		int Nsmpls_rd = 0, Nrows = 0, Ncols = 0; 
		double delta_f; 
		std::vector<double> frq_data; 
		std::vector<std::vector<double>> fft_data; 

		vecut::read_into_vector(freq_file, frq_data, Nsmpls_rd); // read the frequency data back into memory
		delta_f = frq_data[1 + (Nsmpls_rd / 2)]; 
		
		vecut::read_into_matrix(fft_file, fft_data, Nrows, Ncols); // read the computed FFT data back into memory

		std::cout << Nsmpls_rd << " have been read from " << freq_file << "\n"; 
		std::cout << "Frequency sample spacing: "<< delta_f <<"\n";
		std::cout << Nrows << " rows have been read from " << fft_file << "\n";
		std::cout << Ncols << " cols have been read from " << fft_file << "\n\n";

		// store the computed FFT in a single array
		nn = Nsmpls_rd; 
		std::vector<double> fftvals(2 * Nsmpls_rd, 0.0); 

		int count = 0; 
		for (int i = 0; i < Nsmpls_rd; i++) {
			fftvals[count] = fft_data[i][0]; 
			fftvals[count+1] = fft_data[i][1]; 
			count += 2; 
		}

		fft inverse_calc; 

		isign = -1;
		FMT_DATA = false;

		inverse_calc._four1(fftvals, nn, isign, FMT_DATA);

		inverse_calc.output_data(fftvals, delta_f, ift_file, dottxt, isign);

		spctdata.clear(); timedata.clear(); fftvals.clear(); frq_data.clear(); 

	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void testing::real_ft_test()
{
	// test the operation of the real FT code
	// R. Sheehan 6 - 7 - 2022

	try {
		std::string the_dir = "c:\\users\\robertsheehan\\Research\\Notes\\FFT\\Examples";
		useful_funcs::set_directory(the_dir);

		// sine wave example
		int Nsmpls = 4096;

		std::string func_str = "Heaviside";
		std::string timefile = func_str + "_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data
		std::string spctfile = func_str + "_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data

		std::vector<double> timedata; int ntimes = 0;
		std::vector<double> spctdata; int nspct = 0;

		vecut::read_into_vector(timefile, timedata, ntimes);

		vecut::read_into_vector(spctfile, spctdata, nspct);

		double delta_t = timedata[1] - timedata[0];

		unsigned long nn = nspct;

		int isign = +1;

		fft calc;

		calc._realft(spctdata, nn, isign); 

		calc.output_pos_data(spctdata, delta_t, spctfile, dottxt);

		spctdata.clear(); timedata.clear(); 
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void testing::two_ft_test()
{
	// test the operation of the two FT code
	// R. Sheehan 6 - 7 - 2022

	try {
		std::string the_dir = "c:\\users\\robertsheehan\\Research\\Notes\\FFT\\Examples";
		useful_funcs::set_directory(the_dir);

		// sine wave example
		int Nsmpls = 4096;

		std::string func_str = "Sine_Wave";
		std::string timefile1 = func_str + "_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data
		std::string spctfile1 = func_str + "_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data

		std::vector<double> timedata1; int ntimes1 = 0;
		std::vector<double> spctdata1; int nspct1 = 0;

		vecut::read_into_vector(timefile1, timedata1, ntimes1);

		vecut::read_into_vector(spctfile1, spctdata1, nspct1);

		func_str = "Exponential";
		std::string timefile2 = func_str + "_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data
		std::string spctfile2 = func_str + "_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data

		std::vector<double> timedata2; int ntimes2 = 0;
		std::vector<double> spctdata2; int nspct2 = 0;

		vecut::read_into_vector(timefile2, timedata2, ntimes2);

		vecut::read_into_vector(spctfile2, spctdata2, nspct2);

		double delta_t = timedata1[1] - timedata1[0];

		unsigned long nn = nspct1;

		int isign = +1;

		std::vector<double> fft1(2*nspct1, 0.0); 

		std::vector<double> fft2(2*nspct2, 0.0); 

		fft calc; 

		calc._twofft(spctdata1, spctdata2, fft1, fft2, nspct1); 

		//calc.output_data(fft1, delta_t, spctfile1, dottxt); 

		//calc.output_data(fft2, delta_t, spctfile2, dottxt); 

		spctdata1.clear(); spctdata2.clear(); fft1.clear(); fft2.clear(); timedata1.clear(); timedata2.clear(); 
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void testing::convolution_test()
{
	// test the convolution calculation function
	// R. Sheehan 23 - 8 - 2022

	std::string the_dir = "c:\\users\\robertsheehan\\Research\\Notes\\FFT\\Convolutions";
	useful_funcs::set_directory(the_dir);

	// sine wave example
	int Nsmpls = 4096;

	double Lt = 10, x0 = 5, w = 1; 

	tophat_alt(Nsmpls, Lt, x0, w); 

	std::vector<double> timedata; int ntimes = 0;
	std::vector<double> spctdata; int nspct = 0;

	std::string func_str = "Tophat_Alt";
	std::string timefile = func_str + "_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data
	std::string spctfile = func_str + "_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data

	vecut::read_into_vector(timefile, timedata, ntimes);

	vecut::read_into_vector(spctfile, spctdata, nspct);

	double delta_t = timedata[1] - timedata[0];

	unsigned long nn = nspct;

	int isign = +1;

	std::vector<double> ans(nn, 0.0); 
	std::vector<double> respns(spctdata); 

	convol_deconvol calc; 

	calc._convlv(spctdata, nn, respns, nn, isign, ans);

	std::string conv_file = func_str + "_Convolution_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for Convolution data

	std::ofstream write(conv_file, std::ios_base::out, std::ios_base::trunc);

	if (write.is_open()) {

		vecut::wrap_around_conversion(ans); 

		for (size_t i = 0; i < ans.size(); i ++) {
			write << std::setprecision(10) << ans[i] << "\n";
		}

		write.close(); 
	}
}

void testing::sine_wave(int Nsmpls, double Lt, double f1, double f2)
{
	// compute the trace of a multi-frequency sine-wave over fixed time and no. samples
	// R. Sheehan 4 - 7 - 2022

	try {
	
		bool c1 = Nsmpls > 0 ? true : false; 
		bool c2 = Lt > 0 ? true : false; 
		bool c3 = f1 > 0 || f2 > 0 ? true : false; 
		bool c10 = c1 && c2 && c3;

		if (c10) {
			std::string time_file = "Sine_Wave_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data

			std::string data_file = "Sine_Wave_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data

			std::ofstream write1(time_file, std::ios_base::out, std::ios_base::trunc);

			std::ofstream write2(data_file, std::ios_base::out, std::ios_base::trunc);

			if (write1.is_open() && write2.is_open()) {

				double t0 = 0, delta_t = Lt / Nsmpls;
				double omega1 = f1 > 0.0 ? Two_PI * f1 : 0.0;
				double omega2 = f2 > 0.0 ? Two_PI * f2 : 0.0;

				for (int i = 0; i < Nsmpls; i++) {
					write1 << std::setprecision(10) << t0 << "\n";
					write2 << std::setprecision(10) << sin(omega1 * t0) + sin(omega2 * t0) << "\n";
					t0 += delta_t; 
				}

				write1.close();
				write2.close();
			}
		}
		else {
			std::string reason;
			reason = "Error: void testing::sine_wave(int Nsmpls, double Lt, double f1, double f2)\n";
			
			throw std::invalid_argument(reason);
		}

	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void testing::cosine_wave(int Nsmpls, double Lt, double f1, double f2)
{
	// compute the trace of a multi-frequency cosine-wave over fixed time and no. samples
	// R. Sheehan 4 - 7 - 2022

	try {

		bool c1 = Nsmpls > 0 ? true : false;
		bool c2 = Lt > 0 ? true : false;
		bool c3 = f1 > 0 || f2 > 0 ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			std::string time_file = "Cosine_Wave_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data

			std::string data_file = "Cosine_Wave_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data

			std::ofstream write1(time_file, std::ios_base::out, std::ios_base::trunc);

			std::ofstream write2(data_file, std::ios_base::out, std::ios_base::trunc);

			if (write1.is_open() && write2.is_open()) {

				double t0 = 0, delta_t = Lt / Nsmpls;
				double omega1 = f1 > 0.0 ? Two_PI * f1 : 0.0;
				double omega2 = f2 > 0.0 ? Two_PI * f2 : 0.0;

				for (int i = 0; i < Nsmpls; i++) {
					write1 << std::setprecision(10) << t0 << "\n";
					write2 << std::setprecision(10) << cos(omega1 * t0) + cos(omega2 * t0) << "\n";
					t0 += delta_t;
				}

				write1.close();
				write2.close();
			}
		}
		else {
			std::string reason;
			reason = "Error: void testing::cosine_wave(int Nsmpls, double Lt, double f1, double f2)\n";

			throw std::invalid_argument(reason);
		}

	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void testing::gaussian(int Nsmpls, double Lt, double mean, double stdev)
{
	// compute the trace of a Gaussian over fixed time and no. samples
	// Gauss(x) = A exp( -( x - mean )^{2} / ( 2 stdev^{2} ) )
	// R. Sheehan 4 - 7 - 2022

	try {

		bool c1 = Nsmpls > 0 ? true : false;
		bool c2 = Lt > 0 ? true : false;
		bool c3 = mean > 0 && stdev > 0 ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			std::string time_file = "Gaussian_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data

			std::string data_file = "Gaussian_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for sine-wave data

			std::ofstream write1(time_file, std::ios_base::out, std::ios_base::trunc);

			std::ofstream write2(data_file, std::ios_base::out, std::ios_base::trunc);

			if (write1.is_open() && write2.is_open()) {

				double t0 = 0, delta_t = Lt / Nsmpls;
				double numer = 2.0*template_funcs::DSQR(stdev);
				
				for (int i = 0; i < Nsmpls; i++) {
					double arg = (-1.0*template_funcs::DSQR(t0 - mean))/ numer; // -( x - mean )^{2} / ( 2 stdev^{2} )
					write1 << std::setprecision(10) << t0 << "\n";
					write2 << std::setprecision(10) << exp(arg) << "\n";
					t0 += delta_t;
				}

				write1.close(); write2.close(); 
			}
		}
		else {
			std::string reason;
			reason = "Error: void testing::gaussian(int Nsmpls, double Lt, double mean, double stdev)\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void testing::exponential(int Nsmpls, double Lt, double centre, double decay)
{
	// compute the trace of a decaying exponential over fixed time and no. samples
	// Exp(x) = A exp( -decay | x - centre | )
	// R. Sheehan 4 - 7 - 2022

	try {

		bool c1 = Nsmpls > 0 ? true : false;
		bool c2 = Lt > 0 ? true : false;
		bool c3 = decay > 0 && centre > 0 ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			std::string time_file = "Exponential_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for exponential-wave data

			std::string data_file = "Exponential_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for exponential-wave data

			std::ofstream write1(time_file, std::ios_base::out, std::ios_base::trunc);

			std::ofstream write2(data_file, std::ios_base::out, std::ios_base::trunc);

			if (write1.is_open() && write2.is_open()) {

				double t0 = 0, delta_t = Lt / Nsmpls;
				
				for (int i = 0; i < Nsmpls; i++) {
					double arg = -decay * abs(t0 - centre); // -decay | x - centre |
					write1 << std::setprecision(10) << t0 << "\n";
					write2 << std::setprecision(10) << exp(arg) << "\n";
					t0 += delta_t;
				}

				write1.close(); write2.close(); 
			}
		}
		else {
			std::string reason;
			reason = "Error: void testing::exponential(int Nsmpls, double Lt, double centre, double decay)\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void testing::tophat(int Nsmpls, double Lt, double centre, double width)
{
	// compute the trace of a top-hat function over fixed time and no. samples
	// TH(x) = 1 for |x-c| < width, -1 ow
	// R. Sheehan 4 - 7 - 2022

	try {

		bool c1 = Nsmpls > 0 ? true : false;
		bool c2 = Lt > 0 ? true : false;
		bool c3 = width > 0 && centre > 0 ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			std::string time_file = "Tophat_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for exponential-wave data

			std::string data_file = "Tophat_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for exponential-wave data

			std::ofstream write1(time_file, std::ios_base::out, std::ios_base::trunc);

			std::ofstream write2(data_file, std::ios_base::out, std::ios_base::trunc);

			if (write1.is_open() && write2.is_open()) {

				double t0 = 0, delta_t = Lt / Nsmpls;

				for (int i = 0; i < Nsmpls; i++) {
					double val = abs(t0-centre) <= width ? 1.0 : -1.0;
					write1 << std::setprecision(10) << t0 << "\n";
					write2 << std::setprecision(10) << val << "\n";
					t0 += delta_t;
				}

				write1.close(); write2.close();
			}
		}
		else {
			std::string reason;
			reason = "Error: void testing::tophat(int Nsmpls, double Lt, double centre, double width)\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void testing::tophat_alt(int Nsmpls, double Lt, double centre, double width)
{
	// compute the trace of a top-hat function over fixed time and no. samples
	// TH(x) = 1 for |x-c| < width, 0 ow
	// R. Sheehan 4 - 7 - 2022

	try {

		bool c1 = Nsmpls > 0 ? true : false;
		bool c2 = Lt > 0 ? true : false;
		bool c3 = width > 0 && centre > 0 ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			std::string time_file = "Tophat_Alt_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for tophat-wave data

			std::string data_file = "Tophat_Alt_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for tophat-wave data

			std::ofstream write1(time_file, std::ios_base::out, std::ios_base::trunc);

			std::ofstream write2(data_file, std::ios_base::out, std::ios_base::trunc);

			if (write1.is_open() && write2.is_open()) {

				double t0 = 0, delta_t = Lt / Nsmpls;

				for (int i = 0; i < Nsmpls; i++) {
					double val = abs(t0 - centre) <= width ? 1.0 : 0.0;
					write1 << std::setprecision(10) << t0 << "\n";
					write2 << std::setprecision(10) << val << "\n";
					t0 += delta_t;
				}

				write1.close(); write2.close();
			}
		}
		else {
			std::string reason;
			reason = "Error: void testing::tophat(int Nsmpls, double Lt, double centre, double width)\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void testing::heaviside(int Nsmpls, double Lt, double centre)
{
	// compute the trace of a Heaviside function over fixed time and no. samples
	// H(x) = 1 for |x-c| < c, 0 ow
	// R. Sheehan 4 - 7 - 2022

	try {

		bool c1 = Nsmpls > 0 ? true : false;
		bool c2 = Lt > 0 ? true : false;
		bool c3 = centre > 0 ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			std::string time_file = "Heaviside_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for exponential-wave data

			std::string data_file = "Heaviside_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for exponential-wave data

			std::ofstream write1(time_file, std::ios_base::out, std::ios_base::trunc);

			std::ofstream write2(data_file, std::ios_base::out, std::ios_base::trunc);

			if (write1.is_open() && write2.is_open()) {

				double t0 = 0, delta_t = Lt / Nsmpls;

				for (int i = 0; i < Nsmpls; i++) {
					double val = t0 >= centre ? 1.0 : 0;
					write1 << std::setprecision(10) << t0 << "\n";
					write2 << std::setprecision(10) << val << "\n";
					t0 += delta_t;
				}

				write1.close(); write2.close();
			}
		}
		else {
			std::string reason;
			reason = "Error: void testing::heaviside(int Nsmpls, double Lt, double centre)\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void testing::signum(int Nsmpls, double Lt, double centre)
{
	// compute the trace of a signum function over fixed time and no. samples
	// H(x) = 1 for |x-c| < c, -1 ow
	// R. Sheehan 4 - 7 - 2022

	try {

		bool c1 = Nsmpls > 0 ? true : false;
		bool c2 = Lt > 0 ? true : false;
		bool c3 = centre > 0 ? true : false;
		bool c10 = c1 && c2 && c3;

		if (c10) {
			std::string time_file = "Signum_Time_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for exponential-wave data

			std::string data_file = "Signum_Data_Nsmpls_" + template_funcs::toString(Nsmpls) + dottxt; // filename for exponential-wave data

			std::ofstream write1(time_file, std::ios_base::out, std::ios_base::trunc);

			std::ofstream write2(data_file, std::ios_base::out, std::ios_base::trunc);

			if (write1.is_open() && write2.is_open()) {

				double t0 = 0, delta_t = Lt / Nsmpls;

				for (int i = 0; i < Nsmpls; i++) {
					double val = t0 >= centre ? 1.0 : -1.0;
					write1 << std::setprecision(10) << t0 << "\n";
					write2 << std::setprecision(10) << val << "\n";
					t0 += delta_t;
				}

				write1.close(); write2.close();
			}
		}
		else {
			std::string reason;
			reason = "Error: void testing::signum(int Nsmpls, double Lt, double centre)\n";

			throw std::invalid_argument(reason);
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}

void testing::lineshapes()
{
	// compute the FFT of various measured lineshapes
	// R. Sheehan 23 - 2 - 2024

	try {
		std::string the_dir = "C:\\Users\\robertsheehan\\Research\\Laser_Physics\\Linewidth\\Data\\LCR_DSHI_NKT_T_35_D_50\\Beat_Note_Lineshapes\\";
		useful_funcs::set_directory(the_dir);

		// generate the vector for fbeat vals
		int loop_length = 10, Nbeats = 9, delta_f = 80, f_val = delta_f, n_cols, n_rows;
		unsigned long N_spctr_data;
		double frq_spac; 
		std::string filename; 
		for (int i = 0; i < Nbeats; i++) {
			filename = "Lineshape_I_200_D_" + template_funcs::toString(loop_length) + "_fb_" + template_funcs::toString(f_val) + "_fspan_150" + dottxt;
			
			if (useful_funcs::file_exists(filename)) {
				std::cout << "Processing file: " << filename << "\n";

				n_cols = 0;
				n_rows = 0;
				std::vector<std::vector<double>> the_data;

				vecut::read_into_matrix(filename, the_data, n_rows, n_cols);

				//std::cout << "Array Size\n";
				//std::cout << "n_rows: " << n_rows << " , n_cols: " << n_cols << "\n";

				// perform the FFT calculation and output the FFT data
				fft calc;

				std::vector<double> spctr_data(the_data[1]);
				frq_spac = (the_data[0][1] - the_data[0][0]); // frequency spacing in the original data set
				N_spctr_data = n_cols; // no. data points in the original data set

				calc._four1(spctr_data, N_spctr_data);

				calc.output_data(spctr_data, frq_spac, filename, dottxt);

				the_data.clear(); spctr_data.clear();
			}
			else {
				std::cout << "Cannot locate file: " << filename << "\n"; 
			}

			f_val += delta_f;
		}
	}
	catch (std::invalid_argument& e) {
		useful_funcs::exit_failure_output(e.what());
		exit(EXIT_FAILURE);
	}
}