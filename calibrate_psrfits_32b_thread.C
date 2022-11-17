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

    //double StokesVect[4];
    Vector<4, double> StokesVect;
    Vector<4, double> CalStokesVect;
    Matrix<4,4,double> response;
    int nchan = args->nchan;

    std::complex <double> L, arg;
    
    for (int isamp=0; isamp<args->totnpts; isamp++) {
        iptr = (uint8_t *) &args->pfiraw[isamp*npol*nchan + args->ithread*args->nchan_per_thread];

	for (int ichan=0; ichan<args->nchan_per_thread; ichan++) {
	    // Check if zapped channel, then skip
	    //if (args->weights[args->ithread*args->nchan+ichan] == 0){
	    //    iptr++;
	    //continue;
	    //}
	    
	    optr = (float *) &args->pfraw[(isamp*npol*nchan + args->ithread*args->nchan_per_thread + ichan)* 4];
	    
	    response = args->response[ichan + args->ithread*args->nchan_per_thread];
	    
	    I = iptr;
	    Q = (int8_t *)iptr + nchan;
	    U = (int8_t *)iptr + 2*nchan;
	    V = (int8_t *)iptr + 3*nchan;
	    *V *= -1;

	    StokesVect[0] = *I; StokesVect[1] = *Q; StokesVect[2] = *U; StokesVect[3] = *V;

	    CalStokesVect = response * StokesVect;

	    L = std::complex<double>(CalStokesVect[1], CalStokesVect[2]);
	    arg = std::complex<double>(0, 2*(args->rcvr_sa - M_PI/4.));

	    CalStokesVect[1] = real(L*arg); CalStokesVect[2] = imag(L*arg);
	    
	    iptr++;

	    *optr = (float) CalStokesVect[0]; optr += nchan;
	    *optr = (float) CalStokesVect[1]; optr += nchan;
	    *optr = (float) CalStokesVect[2]; optr += nchan;
	    *optr = (float) CalStokesVect[3];
	    
	}
	//exit(-1);
    }
    pthread_exit(NULL);

}
