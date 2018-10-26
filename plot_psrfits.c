#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fitsio.h>
#include <pthread.h>
#include "psrfits.h"
#include "plot_psrfits_cmd.h"
#include "plot_psrfits.h"
#include "dedisperse.h"
#include "cpgplot.h"

#if HAVE_PRESTO
#include "mask.h"
#endif

void reorder_data(unsigned char* outbuf, unsigned char *inbuf, int nband, int nspec, int npol, int nchan, int nbits) {
    int band, spec, pol, inoff = 0, outoff = 0;
    int spband = nspec * npol * nchan;
    int spspec = npol * nchan;

    for (spec = 0 ; spec < nspec ; spec++) {
	for (pol = 0 ; pol < npol ; pol++) {
	    for (band = 0 ; band < nband ; band++) {
		inoff = (band * spband + pol * nchan + spec * spspec) * nbits/8.;
		memcpy(outbuf + outoff, inbuf + inoff, nchan * nbits/8.);
		outoff += nchan * nbits/8.;
	    }
	}
    }
}

static void print_percent_complete(int current, int number, int reset)
{
   static int newper = 0, oldper = -1;

   if (reset) {
      oldper = -1;
      newper = 0;
   } else {
      newper = (int) (current / (float) (number) * 100.0);
      if (newper < 0)
         newper = 0;
      if (newper > 100)
         newper = 100;
      if (newper > oldper) {
         printf("\r%3d%% ", newper);
	 if (newper==100) printf("\n");
         fflush(stdout);
         oldper = newper;
      }
   }
}



int main(int argc, char *argv[]) {

    int i, j, k, status=0;
    Cmdline *cmd;

    // Call usage() if we have no command line arguments
    if (argc == 1) {
        Program = argv[0];
        usage();
        exit(0);
    }
    
    // Parse the command line using the excellent program Clig
    cmd = parseCmdline(argc, argv);

    int numfiles = cmd->argc;
    int fbin = cmd->fbin;
    int tbin = cmd->tbin;
    pthread_t threads[numfiles];  // Thread ids
    thread_args fargs[numfiles];  // Arguments passed to the threads

    // -- Read headers of the PSRFITS files --
    // TODO :  Need to tune basefilename depending on the original name file
    // ie, for Nancay, there is an extra _freq_ term in the fits filename : 'nuppi_55529_0355+54_{scan}_{freq}_0001.fits'
    char *adr, filename[128];
    struct psrfits pf;
    for (i=0; i<numfiles; i++) {
	sprintf(filename, "%s", cmd->argv[i]);
	adr = (char *) strrchr(filename, '_');
	sprintf(adr, "\0");
        sprintf(fargs[i].pf.basefilename, "%s", filename);
	fargs[i].pf.filenum = 1;
        fargs[i].pf.filename[0] = '\0';  
	fargs[i].pf.status = 0;
        status = psrfits_open(&fargs[i].pf);
        if (status) fits_report_error(stderr, status);
    }

    // -- Determine lowest and highest frequency, then set up the middle freq --
    // -- Also determine frequency order -- 
    double flow, fhigh, fctr;
    flow = fhigh = fargs[0].pf.hdr.fctr;
    for (i=1; i<numfiles; i++) {
      if (flow > fargs[i].pf.hdr.fctr) flow = fargs[i].pf.hdr.fctr;
      if (fhigh < fargs[i].pf.hdr.fctr) fhigh = fargs[i].pf.hdr.fctr;
    }
    fctr = flow + (fhigh - flow)/2.;
    if (cmd->verboseP) printf("f_ctr_low = %lf   f_ctr_high = %lf   fctr = %lf\n", flow, fhigh, fctr);

    // -- Identify the channels --
    for (i=0; i<numfiles; i++) {
        fargs[i].chan_id = (int)((fargs[i].pf.hdr.fctr - flow)/fargs[i].pf.hdr.BW);
        if (cmd->verboseP) printf("File #%d  : chan_id = %d\n", i, fargs[i].chan_id);
    }

    // -- Check if we are missing subbands --
    int nbands;
    nbands = (int)((fhigh - flow)/fargs[0].pf.hdr.BW) + 1;
    if (nbands != numfiles) printf("Are we missing one or more subband file ? (%d subbands and %d numfiles)\n", nbands, numfiles);
    printf("Nbands = %d\n", nbands);
    

    // -- Check if same source --
    char source[16];
    sprintf(source, "%s", fargs[0].pf.hdr.source);
    if (cmd->verboseP) printf("Source = %s\n", source);
    for (i=1; i<numfiles; i++) {
        if (strcmp(source, fargs[0].pf.hdr.source)!=0) {
	    fprintf(stderr, " File %s has not the same source as file %s (%s and %s)\n", fargs[i].pf.basefilename, fargs[0].pf.basefilename, fargs[i].pf.hdr.source, fargs[0].pf.hdr.source);
	    exit(-1);
	}    
    }

    // -- Check if same number of channels --
    int nchan;
    nchan = fargs[0].pf.hdr.nchan;
    if(cmd->verboseP) printf("Nchan = %d\n", nchan);
    for (i=1; i<numfiles; i++) {
        if (nchan != fargs[i].pf.hdr.nchan) {
	    fprintf(stderr, " File %s has not the same number of channels as file %s (%d and %d)\n", fargs[i].pf.basefilename, fargs[0].pf.basefilename, fargs[i].pf.hdr.nchan, fargs[0].pf.hdr.nchan);
	    exit(-1);
	}    
    }

    // -- Vars --
    int bytes_per_subint, npol;
    double dt = fargs[0].pf.hdr.dt;
    bytes_per_subint = fargs[0].pf.sub.bytes_per_subint;
    npol = fargs[0].pf.hdr.npol;


    sprintf(pf.hdr.obs_mode, "SEARCH");
    pf.tot_rows = pf.N = pf.T = pf.status = 0;
    pf.filenum = 0;
    pf.filename[0] = '\0';
    pf.hdr.orig_nchan = nchan * nbands;
    pf.hdr.nchan = nchan * nbands;
    int nchans = pf.hdr.nchan;
    pf.hdr.fctr = fctr;
    pf.hdr.BW = fargs[0].pf.hdr.BW * nbands;
    int nsblk = fargs[0].pf.hdr.nsblk;
    int nbits = fargs[0].pf.hdr.nbits;
    double df = fargs[0].pf.hdr.df;
    pf.sub.bytes_per_subint = bytes_per_subint * nbands;
    unsigned long long filelen = (PSRFITS_MAXFILELEN_SEARCH<<30);  
    if (cmd->verboseP) printf("bytes per subint = %d filelen = %lld freq_fact = %d nbits=%d\n", pf.sub.bytes_per_subint, filelen, pf.hdr.ds_freq_fact, nbits);
    pf.rows_per_file = filelen / pf.sub.bytes_per_subint;

    // -- For in-memory transpose of data --
    unsigned char *tmpbuf = NULL;
    unsigned char *tmpbuf32 = NULL;
    tmpbuf = (unsigned char *)malloc(pf.sub.bytes_per_subint * sizeof(unsigned char));
    tmpbuf32 = (unsigned char *)malloc(pf.sub.bytes_per_subint * sizeof(float));

    // -- Alloc data buffers for the output PSRFITS files --
    pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    pf.sub.dat_offsets = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
    pf.sub.dat_scales  = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
    pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint * sizeof(unsigned char));
    pf.sub.data = (unsigned char *)malloc(pf.sub.bytes_per_subint * sizeof(float));

    // -- Given the frequency, set pointers to the correct position --
    for (i=0; i<numfiles; i++) {
        fargs[i].pf.sub.dat_freqs = (float *) &pf.sub.dat_freqs[fargs[i].chan_id * nchan];
        fargs[i].pf.sub.dat_weights = (float *) &pf.sub.dat_weights[fargs[i].chan_id * nchan];
        fargs[i].pf.sub.dat_offsets = (float *) &pf.sub.dat_offsets[fargs[i].chan_id * nchan * npol];
        fargs[i].pf.sub.dat_scales  = (float *) &pf.sub.dat_scales[fargs[i].chan_id * nchan * npol];
        fargs[i].pf.sub.rawdata = (unsigned char *) &tmpbuf[fargs[i].chan_id * bytes_per_subint];
        fargs[i].pf.sub.data = (unsigned char *) &tmpbuf32[fargs[i].chan_id * bytes_per_subint * sizeof(float)];
    }

    //
    int cur_sample = -1;
    int bin_cnt = tbin;
    int start_pt=(int)(cmd->start / dt);
    int end_pt=(int)(cmd->end / dt);
    int nbpts = end_pt - start_pt;
    int nbpts_to_plot = nbpts/tbin;
    nbpts += (int) ((delay_from_dm(cmd->dm, fargs[0].pf.hdr.fctr-fargs[0].pf.hdr.BW) - delay_from_dm(cmd->dm, fargs[numfiles-1].pf.hdr.fctr+fargs[0].pf.hdr.BW) ) / (dt));
    int nbpts_remaining = nbpts;

    unsigned char *f;
    float *data_ptr, *data;
    size_t data_size = (nbpts/tbin+1) * nchans/fbin;
    data = (float *)malloc(data_size * sizeof(float));
    memset(data, 0, data_size * sizeof(float));
    data_ptr = &data[0];

#if HAVE_PRESTO
    int  good_padvals = 0;
    float *padvals;
    mask obsmask;
    if (cmd->maskfile) {
    read_mask(cmd->maskfile, &obsmask);
    printf("Read mask information from '%s'\n", cmd->maskfile);
    //for (i=0;i<1024;i++) printf("%d %d\n", i, obsmask.chans[0][i]);
    good_padvals = determine_padvals(cmd->maskfile, &obsmask, &padvals);
    printf("Will zap %d chans and %d subints\n", obsmask.num_zap_chans, obsmask.num_zap_ints);
    }
#endif

    printf("Will malloc %f MB\n", data_size/1024./1024.);

    printf("Will plot samples %d to %d (%d samples) %d samples with dispersion\n", start_pt, end_pt, nbpts_to_plot, nbpts/tbin);

    // -- Loop through the data --
    int cur_row = 0;
    int statsum = 0;
    int ffact;
    float value = 0.0;
    while (!statsum) {

        // -- Create reading threads -- 
	for (i=0; i<numfiles; i++) {
	    status = pthread_create(&threads[i], NULL, &plot_psrfits_thread, &fargs[i]);
	}    

	// -- Waiting for threads to finish --
	//printf("Lauching threads...\n");fflush(stdout);
	for (i=0; i<numfiles; i++) {
	    pthread_join(threads[i], NULL);
		//if (cmd->verboseP) {
		//    for(j=0;j<fargs[i].pf.hdr.nchan; j++) printf("File #%d  Freq[%03d]=%f\n", i, j, fargs[i].pf.sub.dat_freqs[j]);
		//}
	}
 	for (i=0; i<numfiles; i++)  statsum += fargs[i].status;
	if (statsum) break;

	//copy_subint_params(&pf, &fargs[0].pf);
	//printf("nbands = %d nsblk = %d npol = %d nchan = %d nbits = %d\n", nbands, nsblk, npol, nchan, nbits);fflush(stdout);
        //for (i=0; i<nchan/2; i++) { 
	//    printf("Chan %d %f\n", 2*i+128, (float)( (tmpbuf[npol*nchan/2*nsblk + i] & 240) >> 4) );
        //    printf("Chan %d %f\n", 2*i+128+1, (float)((tmpbuf[npol*nchan/2*nsblk + i] & 15)));
	//}     
	//exit(0);
	reorder_data(pf.sub.rawdata, tmpbuf, nbands, nsblk, npol, nchan, nbits);

	for (j=0; j < nsblk; j++) {
	    cur_sample += 1;
	    //printf("cur sample = %d, start_pt = %d\n", cur_sample, start_pt);
	    if (cur_sample < start_pt )
		continue;
	    if (nbpts_remaining) {
	        //printf("cur sample = %d\n", cur_sample);
		f = (unsigned char *) &pf.sub.rawdata[j * nchans/2 * npol];
		ffact = 1;
		for(i=0;i<nchans/2;i++) {
		    value = (float)( (f[i] & 240) >> 4);
#if HAVE_PRESTO		    
		if (cmd->maskfile) {
		    for(k=0; k<obsmask.num_chans_per_int[cur_row]; k++) {
		        //printf("cur_row = %d k=%d cur_chan=%d\n", cur_row, k, 2*i);
		        if (obsmask.chans[cur_row][k] == 2*i) {
			    value = padvals[2*i];
			    //printf("Zap chan %d with val=%f instead of %f\n", 2*i, value, (float)( (f[i] & 240) >> 4));
			    break;
			}    
		    }	
		}	
#endif		    
		    *data_ptr += value;

		    if (fbin==1) data_ptr++;
		    value = (float)(f[i] & 15);

#if HAVE_PRESTO		
		if (cmd->maskfile) {
		    for(k=0; k<obsmask.num_chans_per_int[cur_row]; k++) {
		        if (obsmask.chans[cur_row][k] == 2*i+1) {
			    value = padvals[2*i+1];
			    //printf("Zap chan %d with val=%f\n", 2*i+1, value);
			    break;
			}
		    }
		}
#endif		    
		    *data_ptr += value;

		    if (fbin==1 || ffact == fbin/2 ) {data_ptr++; ffact = 1;}
		    else ffact ++;
		}
		nbpts_remaining--;

		// Binning counter
		bin_cnt--;
		if(bin_cnt) {
		    data_ptr -= nchans / fbin;
		}
		else bin_cnt=tbin;

	    }  
	    //cur_sample += 1;
	    //

	}
	//printf("subint chan %d %d\n", cur_row, obsmask.chans[cur_row][0]);
	cur_row ++;

	//for(i=0;i<nchans;i++) {


	print_percent_complete(fargs[0].pf.tot_rows, fargs[0].pf.rows_per_file, fargs[0].pf.tot_rows == 1 ? 1:0);
    } 

  //for (i=0; i<nchans; i++) printf("Chan %d %f\n", i, data[i]);

    // -- Closing all files --
    for (i=0; i<numfiles; i++) {
        status = psrfits_close(&fargs[i].pf);
    }	

    // X coordinates for the dedispersed time serie
    float *xdata;
    xdata=(float *)malloc(nbpts_to_plot * sizeof(float));
    for(i=0;i<nbpts_to_plot;i++) xdata[i] = cmd->start+i*dt*tbin;


  // -- Normalisation between 0 and 1 for channels --
  float xmin,xmax,ymin,ymax,y_diff,val_min,val_max;
  val_min = data[0];
  val_max = data[0];
  for(i=0;i<nchans/fbin;i++) {
    for(j=0; j<nbpts/tbin; j++) {
      //printf("%u \n", (unsigned int)data[i+j*nchans/fbin]);
      if(data[i+j*nchans/fbin]>val_max) val_max=data[i+j*nchans/fbin];
      if(data[i+j*nchans/fbin]<val_min) val_min=data[i+j*nchans/fbin];
      //if(i==0) printf("%d %u\n", j, (unsigned int) data[i+j*nchans/fbin] );
    }
  }
  printf("Max = %f Min =%f\n", val_max, val_min);

  for(i=0;i<nchans/fbin;i++) {
    if(val_min != val_max) {
      for(j=0;j<nbpts/tbin;j++) 
	data[i+j*nchans/fbin]=(data[i+j*nchans/fbin]-val_min)/(val_max-val_min);      
    } else {
      for(j=0;j<nbpts/tbin;j++) data[i+j*nchans/fbin]=0.0;
    }

  }

  // -- Open PgPlot output
  cpgopen("plot.ps/ps");

  // -- Set window --
  xmin=cmd->start-0.0125*(cmd->end-cmd->start);xmax=cmd->end+0.0125*(cmd->end-cmd->start);
  ymin=pf.sub.dat_freqs[0]; 
  ymax=pf.sub.dat_freqs[nchans-1]; 
  y_diff=ymax-ymin;
  ymin=ymin-0.025*y_diff;ymax=ymax+0.025*y_diff;

  cpgsvp(0.1,0.9,0.1,0.7);
  printf("xmin=%f xmax=%f ymin=%f ymax=%f\n",xmin,xmax,ymin,ymax);
  cpgswin(xmin,xmax,ymin,ymax);

  // -- Matrix --
  float tr[6];
  tr[0]=cmd->start; tr[1]=0.0; tr[2]=dt*tbin;
  tr[3]=pf.hdr.fctr - nchans/2 * df; tr[4]=df*fbin; tr[5]=0.0;
  for (i=0;i<6;i++) printf("tr[%d] = %f\n", i, tr[i]);

  val_min=data[0];val_max=data[0];
  for(i=0;i<nchans/fbin*nbpts/tbin;i++) {
    //printf("%d %u\n", i, (unsigned int) data[i] );
    if (data[i]>val_max) val_max=data[i];
    if (data[i]<val_min) val_min=data[i];
  }

  // -- Dedisperse the data --
  // -- First, construct delays with respect to the highest frequency
  int *idelays = (int *)malloc(nchans/fbin * sizeof(int));
  double freq = 0.0;
  double hifreq = 0.0;

  for (i=nchans-1; i>nchans-fbin-1; i--) hifreq += pf.sub.dat_freqs[i];
  hifreq /= fbin;

  for(i=0;i<nchans/fbin;i++) {
      freq = 0.0;
      for (j=0; j<fbin; j++) freq += pf.sub.dat_freqs[i*fbin+j];
      freq /= fbin;
      idelays[i] = (int) ((delay_from_dm(cmd->dm, freq) - delay_from_dm(cmd->dm, hifreq) ) / (dt*tbin));
      //idelays[nchans/fbin-i-1] = (int) ((delay_from_dm(cmd->dm, freq) - delay_from_dm(cmd->dm, hifreq) ) / (dt*tbin));
      printf("freq = %lf  hifreq = %lf  dt=%lf dm=%f   freq delay (in bins) = %d\n", freq, hifreq, dt,  cmd->dm, idelays[i]);
  }

  // -- Then do the dedispersion
  float *dedis_data = (float *)malloc(nbpts/tbin * sizeof(float));
  memset(dedis_data, 0, nbpts/tbin * sizeof(float));
  int max_delays = idelays[0];
  dedisperse(data, dedis_data, idelays, max_delays, nchans/fbin, nbpts/tbin, 0); // max_delays prevent to go beyond the allocated memory at the lowest frequency


  // -- Plot the 2D Map --
  val_max=1.0; val_min=0.0;
  printf("Data: %f %f %f %f %f %f %f %f\n", data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7]);
  cpggray(data, nchans/fbin, nbpts/tbin, 1, nchans/fbin, 1, nbpts_to_plot, val_max, val_min, tr);
  cpglab("","Frequency (MHz)","");
  cpgbox("BCSNT",0.0,0,"BCSNT",0.0,0);
  cpgmtxt("B",3.0,0.5,0.5,"Time (s)");
  

  // -- Plot the dedispersed time serie -- 
  val_max=0.0; val_min=1.0;
  for(i=0;i<nbpts_to_plot;i++) {
    //printf("%d %f\n", i, dedis_data[i]);
    if (dedis_data[i]>val_max) val_max=dedis_data[i];
    if (dedis_data[i]<val_min) val_min=dedis_data[i];
  }
  cpgsvp(0.1,0.9,0.7,0.9);
  cpgswin(xmin,xmax, val_min, val_max);
  cpgbox("BCST",0.0,0,"BCST",0.0,0);
  cpgline(nbpts_to_plot, xdata, dedis_data);


  // -- Print title --
  char string[128];
  sprintf(string,"%s RA=%s DEC=%s", fargs[0].pf.hdr.source, fargs[0].pf.hdr.ra_str, fargs[0].pf.hdr.dec_str);
  cpgmtxt("T",2.0,0.5,0.5,string);

  sprintf(string,"TBin=%d   FBin=%d   TimeResol=%.1lf\\gms   DM=%.1fpc/cm3", tbin, fbin, dt*(float)tbin*1e6, cmd->dm);
  cpgsch(0.7); cpgmtxt("T",0.5,1.0,1.0,string); cpgsch(1.0);

  cpgclos();


    exit(0);
}
