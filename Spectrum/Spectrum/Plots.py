# python script for plotting computed FFT spectra
# R. Sheehan 9 - 8 - 2019

import os
import sys
import glob

import math
import scipy
import numpy as np
import matplotlib.pyplot as plt

# add path to our file
sys.path.append('c:/Users/robertsheehan/Programming/Python/Common/')
sys.path.append('c:/Users/robertsheehan/Programming/Python/Plotting/')

import Common
import Plotting

MOD_NAME_STR = "Plots" # use this in exception handling messages

# method definitions

def general_fft_plot():
	# make a plot of the general FFT that is computed
	# R. Sheehan 9 - 8 - 2019
	
	FUNC_NAME = ".general_fft_plot()" # use this in exception handling messages
	ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME
	
	try:
		os.chdir('c:/Users/robertsheehan/Research/Notes/FFT/Examples/')

		Nsmpls = 4096
		wavename = 'Unequal_Ampl_Sine_Wave'

		time_file = "%(v2)s_Time_Nsmpls_%(v1)d.txt"%{"v2":wavename, "v1":Nsmpls}		
		spct_file = "%(v2)s_Data_Nsmpls_%(v1)d.txt"%{"v2":wavename, "v1":Nsmpls}

		IFT = False; 

		if IFT:		
			frq_file = "%(v2)s_Data_Nsmpls_%(v1)d_IFT_Frq_data.txt"%{"v2":wavename, "v1":Nsmpls}
			fft_file = "%(v2)s_Data_Nsmpls_%(v1)d_IFT_FFT_data.txt"%{"v2":wavename, "v1":Nsmpls}
		else:
			frq_file = "%(v2)s_Data_Nsmpls_%(v1)d_Frq_data.txt"%{"v2":wavename, "v1":Nsmpls}
			fft_file = "%(v2)s_Data_Nsmpls_%(v1)d_FFT_data.txt"%{"v2":wavename, "v1":Nsmpls}
			fft_wrap = "%(v2)s_Data_Nsmpls_%(v1)d_FFT_data_wrap_around.txt"%{"v2":wavename, "v1":Nsmpls}

		if glob.glob(time_file) and glob.glob(spct_file) and glob.glob(frq_file) and glob.glob(fft_file):
			
			# plot the time series
			time_data = np.loadtxt(time_file, unpack = True)
			spct_data = np.loadtxt(spct_file, unpack = True)
			
			args = Plotting.plot_arg_single()
			
			args.loud = False
			args.x_label = 'Time / time-units'
			args.y_label = 'Signal / volt-units'
			args.marker = 'r-'
			args.fig_name = spct_file.replace('.txt','')
			
			Plotting.plot_single_curve(time_data, spct_data, args)
			
			del time_data; del spct_data; 

			# plot the FFT in wrap-around format
			#spct_data = np.loadtxt(fft_wrap, unpack = True)
			#time_data = np.arange(0, len(spct_data), 1)

			#args.loud = False
			#args.x_label = 'Time / time-units'
			#args.y_label = 'Signal / volt-units'
			#args.marker = 'r-'
			#args.fig_name = 'Signal_FFT_Wrap_Around'
			
			#Plotting.plot_single_curve(time_data, spct_data, args)

			#del time_data; del spct_data; 
			
			# plot the computed FFT
			frq_data = np.loadtxt(frq_file, unpack = True)
			spct_data = np.loadtxt(fft_file, delimiter = ',', unpack = True)

			# need to scale plots according to max{Re{FFT}} or max{Im{FFT}}
			max_re = np.max(spct_data[0]); max_im = np.max(spct_data[1]); 
			max_val = max(max_re, max_im)
			#max_val = np.max(spct_data[2])

			hv_data = []; marks = []; labs = []; 
			hv_data.append([frq_data, spct_data[0] / max_val ]); marks.append(Plotting.labs_lins[0]); labs.append('Re{FFT}')
			hv_data.append([frq_data, spct_data[1] / max_val ]); marks.append(Plotting.labs_lins[1]); labs.append('Im{FFT}')
			#hv_data.append([frq_data, spct_data[2] / max_val ] ); marks.append(Plotting.labs_lins[2]); labs.append('Abs{FFT}')
			#hv_data.append([frq_data, spct_data[3]]); marks.append(Plotting.labs_lins[3]); labs.append('Arg{FFT}')

			args = Plotting.plot_arg_multiple()

			args.loud = True
			args.x_label = 'Frequency / 1/time-units'
			args.y_label = 'Signal FFT'
			args.crv_lab_list = labs
			args.mrk_list = marks
			args.plt_range = [-40, 40, -1, 1] if IFT is False else None
			args.fig_name = fft_file.replace('.txt','')		
			Plotting.plot_multiple_curves(hv_data, args)
			
			del frq_data; del spct_data; del hv_data; del labs; del marks; 
			
		else:
			ERR_STATEMENT = ERR_STATEMENT + "\nInput file not found"
			raise Exception		
	except Exception as e:
		print(ERR_STATEMENT)
		print(e)

def Single_plot():
	
	os.chdir('c:/Users/robertsheehan/Research/Notes/FFT/Examples/')

	#spct_file = "TwoFFT_WAO.txt"
	#spct_file = "TwoFFT_Unwrapped_Neg.txt"
	spct_file = "TwoFFT_Separated.txt"
	#spct_file = "TwoFFT_Equality.txt"


	# plot the FFT in wrap-around format
	spct_data = np.loadtxt(spct_file, delimiter = ',', unpack = True)
	#time_data = np.arange(0, len(spct_data[0]), 1)
	time_data = []
	df = 0.1; f0=0.1; 
	for i in range(0, len(spct_data[0]), 1):
		time_data.append(f0)
		f0 = f0 + 0.5*df

	args = Plotting.plot_arg_single()
	args.loud = True
	args.x_label = 'x-axis'
	args.y_label = 'y-axis'
	args.marker = 'r-'
	args.fig_name = 'F'
	
	Plotting.plot_single_curve(time_data, spct_data[0], args)

	args.marker = 'g-'
	args.fig_name = 'G'
	
	Plotting.plot_single_curve(time_data, spct_data[1], args)

	#args.marker = 'c-'
	#args.fig_name = 'Fn+FN-n'	
	#Plotting.plot_single_curve(time_data, spct_data[0] - spct_data[1], args)

	#args = Plotting.plot_arg_multiple()
	#args.loud = True
	#args.x_label = 'x-axis'
	#args.y_label = 'y-axis'
	##args.mrk_list = ['r-','g-','c-','m-']
	##args.crv_lab_list = ['0','1','2','3']
	#args.mrk_list = ['r-']
	#args.crv_lab_list = ['0']
	##args.marker = 'r-'
	##args.fig_name = 'Signal_FFT_Wrap_Around'
	
	#hv_data = []; 
	#hv_data.append([time_data, spct_data[0]])
	##hv_data.append([time_data, spct_data[1]])
	##hv_data.append([time_data, spct_data[2]])
	##hv_data.append([time_data, spct_data[3]])
	
	#Plotting.plot_multiple_curves(hv_data, args)

	del time_data; del spct_data; 

def laser_fft_plot():
	# make a plot of the FFT of a laser optical spectrum
	# deduce the laser cavity length from the FFT
	# R. Sheehan 9 - 8 - 2019
	
	FUNC_NAME = ".laser_fft_plot()" # use this in exception handling messages
	ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME
	
	try:
		DATA_HOME = 'c:/users/robert/Research/FFT_Data/Examples_2019/'
		
		if os.path.isdir(DATA_HOME):
			os.chdir(DATA_HOME)
			print(os.getcwd())
			
			#wl_file = "ed_laser_wl.csv"
			#pow_file = "ed_laser_pow.csv"		
			
			wl_file = "PM_Laser_WL.csv"
			pow_file = "PM_Laser_Pow.csv"
			
			#wl_file = "WL_Meas_Shallow_I_90.txt"
			#pow_file = "Pow_Meas_Shallow_I_90.txt"
			
			#wl_file = "20C_wave.dat"
			#pow_file = "50mA.dat"
			
			# wl_file = "Wavelength.txt"
			# pow_file = "Spectrum_I_150.txt"
			
			#extension = ".txt"
			extension = ".csv"
			#extension = ".dat"
			
			frq_file = pow_file.replace(extension,"") + "_Frq_data" + extension
			fft_file = pow_file.replace(extension,"") + "_Abs_FFT_data" + extension
			
			if glob.glob(wl_file) and glob.glob(pow_file) and glob.glob(frq_file) and glob.glob(fft_file):
				# plot the time series
				wl_data = np.loadtxt(wl_file, unpack = True)
				spct_data = np.loadtxt(pow_file, unpack = True)
				
				args = Plotting.plot_arg_single()
				
				args.loud = True
				args.x_label = 'Wavelength / nm'
				args.y_label = 'Power / dBm'
				args.marker = 'r-'
				args.fig_name = pow_file.replace(extension,"")
				
				Plotting.plot_single_curve(wl_data, spct_data, args)
				
				del wl_data; del spct_data;
				
				# plot the computed FFT
				frq_data = np.loadtxt(frq_file, unpack = True)
				spct_data = np.loadtxt(fft_file, unpack = True)
				
				lambda_1 = 1500.0
				#n_g = 3.2 # estimate of the material group index
				n_g = 1.3 # estimate of the material group index
				for i in range(0, len(frq_data), 1):
					frq_data[i] = ( (lambda_1**2)/( 2.0*n_g ) )*frq_data[i]/1000.0 # divide by 1000 to convert from nm to um
				
				args.loud = True
				args.x_label = 'Cavity Length / um'
				args.y_label = 'Signal FFT'
				args.marker = 'g-'
				args.fig_name = fft_file.replace(extension,"")
				args.plt_range = [0, 8e+3, 0, 30e+3]
				
				Plotting.plot_single_curve(frq_data, spct_data, args)
				
				del frq_data; del spct_data;
			else:
				ERR_STATEMENT = ERR_STATEMENT + "\nCannot locate input files"
				raise Exception			
		else:
			raise EnvironmentError
	except EnvironmentError:
		print(ERR_STATEMENT)
		print("Cannot locate directory:",DATA_HOME)
	except Exception as e:
		print(ERR_STATEMENT)
		print(e)

# method calls
general_fft_plot()

#Single_plot()

#laser_fft_plot()