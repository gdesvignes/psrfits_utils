#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include <pthread.h>
#include <cstdint>
#include "psrfits.h"
#include "calibrate_psrfits.h"
#include <iostream>
#include <complex>
#include "Pauli.h"

void *calibrate_psrfits_32b_thread(void *_args) {

    uint8_t *uptr;
    int status = 0, npol=4;
    uint8_t *I;
    int8_t *Q, *U, *V;

    // -- Get arguments --
    thread_args *args = (thread_args *) _args;

    // Input pointers
    uint8_t *iptr;

    // Output pointers
    float *optr;

    Vector<4, double> StokesVect;
    Vector<4, double> CalStokesVect;
    //Matrix<4,4,double> response;
    
    int nchan = args->nchan;

    std::vector<Matrix<4,4,double>> response;
    std::complex <double> L, arg;

    Matrix<4,4,double> Mzero;

    for (int ichan=args->ithread*args->nchan_per_thread; ichan < (args->ithread+1)*args->nchan_per_thread; ichan++) {
      if (args->is_chan_valid[ichan] == 1)
	response.push_back( inv(args->response[ichan] * args->MPA) );
      else
	response.push_back(Mzero);
    }

    int osamp;
    for (int isamp=0; isamp<args->totnpts; isamp++) {
        iptr = (uint8_t *) &args->pfiraw[isamp*npol*nchan + args->ithread*args->nchan_per_thread];

	osamp = isamp/args->ds_time_fact; // integer division truncates the decimals (if any)
	
	for (int ichan=0; ichan<args->nchan_per_thread; ichan++) {
	    
	    optr = (float *) &args->pfraw[(osamp*npol*nchan + args->ithread*args->nchan_per_thread + ichan)* 4];
	    
	    //response = args->response[ichan + args->ithread*args->nchan_per_thread];
	    
	    I = iptr;
	    Q = (int8_t *)iptr + nchan;
	    U = (int8_t *)iptr + 2*nchan;
	    V = (int8_t *)iptr + 3*nchan;
	    *V *= -1;

	    StokesVect[0] = *I; StokesVect[1] = *Q; StokesVect[2] = *U; StokesVect[3] = *V;

	    CalStokesVect = response[ichan] * StokesVect;

	    L = std::complex<double>(CalStokesVect[1], CalStokesVect[2]);
	    arg = std::complex<double>(0, 2*(args->rcvr_sa - M_PI/4.));

	    CalStokesVect[1] = real(L*arg); CalStokesVect[2] = imag(L*arg);
	    
	    iptr++;

	    *optr += (float) CalStokesVect[0]; optr += nchan;
	    *optr += (float) CalStokesVect[1]; optr += nchan;
	    *optr += (float) CalStokesVect[2]; optr += nchan;
	    *optr += (float) CalStokesVect[3];
	    
	}
    }
    pthread_exit(NULL);

}
