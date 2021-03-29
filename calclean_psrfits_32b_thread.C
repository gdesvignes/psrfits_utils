#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include <pthread.h>
#include "psrfits.h"
#include "calclean_psrfits.h"
//#include "ekstrom.h"

#include <iostream>

void *calclean_psrfits_32b_thread(void *_args) {

    int status = 0;
	float value, median;    

    // -- Get arguments --
    thread_args *args = (thread_args *) _args;
	
	// Input pointers
	float *pfiraw = (float *) args->pfiraw;
	pfiraw += args->ithread*args->nchan_per_thread;
	
	// Output pointers
	float *pfraw = (float *) args->pfraw;
	pfraw += args->ithread*args->nchan_per_thread;
	
	unsigned long index;

	float *median_ar = (float*) malloc(args->nchan_per_thread * args->npol * sizeof(float));

	for (int isamp=0; isamp<args->totnpts; isamp++) {
	  for (int ipol=0; ipol<args->npol; ipol++) {
		for (int ichan=0; ichan<args->nchan_per_thread; ichan++) {

		  value = *pfiraw;
		  if (isamp%args->nskip==0)
			median_ar[ipol * args->nchan_per_thread + ichan] = args->rm[ipol * args->nchan_per_thread + ichan].get_value(value);

		  *pfraw = (value - median_ar[ipol * args->nchan_per_thread + ichan]); // Assumes output is 8bits

		  pfiraw++;
		  pfraw++;

		}
		pfiraw += (args->nchan-args->nchan_per_thread);
		pfraw += (args->nchan-args->nchan_per_thread);
		/*
		  for (int ipol=1; ipol<args->npol; ipol++) {
		  for (int ichan=0; ichan<args->nchan_per_thread; ichan++) {
		  // GD value = (float) *pfiraw_s;
		  //if (isamp%args->nskip==0)
		  //  median_ar[ipol * args->nchan_per_thread + ichan] = args->rm[ipol * args->nchan_per_thread + ichan].get_value(value);
		  // GD *pfraw_s = (value - median_ar[ipol * args->nchan_per_thread + ichan]);
		  *pfraw = *pfiraw;
		  pfiraw++;
		  pfraw++;
		}

		pfiraw += (args->nchan-args->nchan_per_thread);
		pfraw += (args->nchan-args->nchan_per_thread);
		}*/
	  } 
	}
	free(median_ar);

    pthread_exit(NULL);

}

