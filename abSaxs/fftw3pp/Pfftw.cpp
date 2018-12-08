/*
 * Pfftw.cpp
 *
 *  Created on: Dec 7, 2018
 *      Author: marchi
 */

#include "../fftw3pp/Pfftw.h"

namespace fftwpp {
std::ifstream Pfftw::ifWisdom;
std::ofstream Pfftw::ofWisdom;
bool Pfftw::Wise=false;
const double Pfftw::twopi=2.0*acos(-1.0);

// User settings:
unsigned int Pfftw::effort=FFTW_MEASURE;
const char *Pfftw::WisdomName="wisdom3.txt";

} /* namespace fftwpp */
