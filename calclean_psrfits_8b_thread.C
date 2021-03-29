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

void *calclean_psrfits_8b_thread(void *_args) {

    int status = 0;
	float value, median;    

    // -- Get arguments --
    thread_args *args = (thread_args *) _args;
	
	// Input pointers
	unsigned char *pfiraw;
	char *pfiraw_s;
	
	// Output pointers
	unsigned char *pfraw = (unsigned char*) args->pfraw;
	char *pfraw_s = (char *) args->pfraw;
	// if 32bits output, all is float --
	if (args->nbits==32) {
	  float *pfraw = (float *) args->pfraw;
	}
	
	unsigned long index;

	float *median_ar = (float*) malloc(args->nchan_per_thread * args->npol * sizeof(float));
	
	pfiraw = (unsigned char *) &args->pfiraw[args->ithread*args->nchan_per_thread];
	pfiraw_s = (char *) &args->pfiraw[args->ithread*args->nchan_per_thread];

	pfraw += args->ithread*args->nchan_per_thread;
	pfraw_s += args->ithread*args->nchan_per_thread;
	
	for (int isamp=0; isamp<args->totnpts; isamp++) {
	  for (int ichan=0; ichan<args->nchan_per_thread; ichan++) {

		value = (float) *pfiraw;
		if (isamp%args->nskip==0)
		  median_ar[ichan] = args->rm[ichan].get_value(value);

		*pfraw = (value - median_ar[ichan] + 127); // Assumes output is 8bits

		//if (ichan==451) {
		  //printf("isamp=%d ithread=%d totsamp=%d ichan=%d ipol=0 value = %f median = %f %p %p\n", isamp, args->ithread, args->totnpts, ichan, value, median_ar[ichan], pfiraw, pfraw);
		  //printf("%d %d\n", (int)value, (int) *pfraw);
		//}
		
		pfiraw++;
		pfiraw_s++;
		pfraw++;
		pfraw_s++;
	  }
	  pfiraw += (args->nchan-args->nchan_per_thread);
	  pfiraw_s += (args->nchan-args->nchan_per_thread);
	  pfraw += (args->nchan-args->nchan_per_thread);
	  pfraw_s += (args->nchan-args->nchan_per_thread);

	  for (int ipol=1; ipol<args->npol; ipol++) {
		for (int ichan=0; ichan<args->nchan_per_thread; ichan++) {
		  // GD value = (float) *pfiraw_s;
		  //if (isamp%args->nskip==0)
		  //  median_ar[ipol * args->nchan_per_thread + ichan] = args->rm[ipol * args->nchan_per_thread + ichan].get_value(value);
		  // GD *pfraw_s = (value - median_ar[ipol * args->nchan_per_thread + ichan]);
		  *pfraw_s = *pfiraw_s;
		  pfiraw++;
		  pfiraw_s++;
		  pfraw++;
		  pfraw_s++;
		}

		pfiraw += (args->nchan-args->nchan_per_thread);
		pfiraw_s += (args->nchan-args->nchan_per_thread);
		pfraw += (args->nchan-args->nchan_per_thread);
		pfraw_s += (args->nchan-args->nchan_per_thread);
	  }

	  //pfiraw += (args->nchan-args->nchan_per_thread)*args->npol;
	  //pfiraw_s += (args->nchan-args->nchan_per_thread)*args->npol;
	  //pfraw += (args->nchan-args->nchan_per_thread)*args->npol;
	  //pfraw_s += (args->nchan-args->nchan_per_thread)*args->npol;
	  
	}

	free(median_ar);

    pthread_exit(NULL);

}

