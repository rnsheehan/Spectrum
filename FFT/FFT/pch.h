// pch.h: This is a precompiled header file.
// Files listed below are compiled only once, improving build performance for future builds.
// This also affects IntelliSense performance, including code completion and many code browsing features.
// However, files listed here are ALL re-compiled if any one of them is updated between builds.
// Do not add files here that you will be updating frequently as this negates the performance advantage.

#ifndef PCH_H
#define PCH_H

// add headers that you want to pre-compile here
#include "framework.h"

#include <cstdlib>
#include <iostream> // cout, cin, cerr
#include <iomanip> // setw, setprecision, time
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <limits.h>
#include <utility>
#include <direct.h> // needed for _mkdir, _chdir

// Constants
static const double EPS = (3.0e-12);
static const double p = (atan(1.0)); // pi / 4
static const double Two_PI = (8.0 * p); // 2 pi
static const double PI = (4.0 * p); // pi
static const double PI_2 = (2.0 * p); // pi / 2
static const double PI_3 = ((4.0 / 3.0) * p); // pi / 3
static const double PI_4 = (p); // pi / 4
static const double PI_5 = ((4.0 / 5.0) * p); // pi / 5
static const double PI_6 = ((2.0 / 3.0) * p); // pi / 6 

static const int MAX_PATH_LENGTH = 250; // max. length for a directory in Windows OS

static const std::string empty_str = "";
static const std::string dottxt = ".txt";
static const std::string dotcsv = ".csv";

#include "Templates.h"
#include "Useful.h"
#include "Vector_Utils.h"
#include "FFT_ALG.h"
#include "FFT.h"

#endif //PCH_H
