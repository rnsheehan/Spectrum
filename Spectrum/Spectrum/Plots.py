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
sys.path.append('c:/Users/Robert/Programming/Python/Common/')
sys.path.append('c:/Users/Robert/Programming/Python/Plotting/')

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
		time_file = "Time_Data.txt"		
		spct_file = "Spec_Data.txt"
		
		frq_file = "Spec_Data_Frq_data.txt"
		fft_file = "Spec_Data_Abs_FFT_data.txt"
		
		if glob.glob(time_file) and glob.glob(spct_file) and glob.glob(frq_file) and glob.glob(fft_file):
			
			# plot the time series
			time_data = np.loadtxt(time_file, unpack = True)
			spct_data = np.loadtxt(spct_file, unpack = True)
			
			args = Plotting.plot_arg_single()
			
			args.loud = True
			args.x_label = 'Time / s'
			args.y_label = 'Signal'
			args.marker = 'r-'
			args.fig_name = 'Signal_Data'
			
			Plotting.plot_single_curve(time_data, spct_data, args)
			
			del time_data; del spct_data; 
			
			# plot the computed FFT
			frq_data = np.loadtxt(frq_file, unpack = True)
			spct_data = np.loadtxt(fft_file, unpack = True)
			
			args.loud = True
			args.x_label = 'Frequency / Hz'
			args.y_label = 'Signal FFT'
			args.marker = 'g-'
			args.fig_name = 'Signal_FFT'
			
			Plotting.plot_single_curve(frq_data, spct_data, args)
			
			del frq_data; del spct_data; 
			
		else:
			ERR_STATEMENT = ERR_STATEMENT + "\nInput file not found"
			raise Exception		
	except Exception as e:
		print(ERR_STATEMENT)
		print(e)

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
#general_fft_plot()

laser_fft_plot()


	
	