#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include <pthread.h>
#include "psrfits.h"
#include "calclean_psrfits.h"

#include "presto.h"

void *calclean_psrfits_thread(void *_args) {

    int status = 0;
    int numbins;
    

    // -- Get arguments --
    thread_args *args = (thread_args *) _args;

    pthread_mutex_lock(&lock_calclean);


    fcomplex *fft;
    if (args->clean) {
    
	/* FFT it */
	realfft(args->data, args->numdata, -1);
	fft = (fcomplex *) args->data;
	numbins = args->numdata / 2;
	//printf("done.\n");

	/* De-redden it */
	deredden(fft, numbins);
	
	/* FFT back */
	realfft(args->data, args->numdata, 1);

    }


    pthread_mutex_unlock(&lock_calclean);
    pthread_exit(NULL);

}

