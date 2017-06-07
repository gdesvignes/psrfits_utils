#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <getopt.h>
#include "psrfits.h"
#include "ascii_header.h"

#ifndef DEGTORAD
#define DEGTORAD 0.017453292519943295769236907684886127134428718885417
#endif
#ifndef RADTODEG
#define RADTODEG 57.29577951308232087679815481410517033240547246656
#endif

#define DADA_HEADER_SIZE (4096)

void usage() {
  printf(
		 "Usage: fold_psrfits [options] input_filename_base\n"
		 "Options:\n"
		 "  -h, --help               Print this\n"
		 "  -f fff, --freq=ffff      Set the centre frequency (MHz)\n"
		 "  -i, --invert             Invert band\n"
		 );
}

void dec2hms(char *out, double in, int sflag) {
    int sign = 1;
    char *ptr = out;
    int h, m;
    double s;
    if (in<0.0) { sign = -1; in = fabs(in); }
    h = (int)in; in -= (double)h; in *= 60.0;
    m = (int)in; in -= (double)m; in *= 60.0;
    s = in;
    if (sign==1 && sflag) { *ptr='+'; ptr++; }
    else if (sign==-1) { *ptr='-'; ptr++; }
    sprintf(ptr, "%2.2d:%2.2d:%07.4f", h, m, s);
}

int main(int argc, char *argv[]) {

  /* Cmd line */
  static struct option long_opts[] = {
	{"freq",    1, NULL, 'f'},
	{"invert",  0, NULL, 'i'},
	{"help",    0, NULL, 'h'},
	{0,0,0,0}
  };

  int opt, opti;
    int ii, jj;
    double dtmp;
    struct psrfits pf;
    
	
	
    FILE *pfi;
    printf("Filename = %s\n", argv[1]);
    pfi = fopen(argv[1], "r");
	
	
	char header[DADA_HEADER_SIZE];
	
	if (fread (header, 1, DADA_HEADER_SIZE, pfi) != DADA_HEADER_SIZE) {
	  fprintf (stderr, "Error getting header from %s\n", argv[1]);
	  exit(-1);
	}
	
	double input_freq;
	bool have_freq=false, have_invert=false;
	while ((opt=getopt_long(argc,argv,"f:ih", long_opts,&opti))!=-1) {
	  switch (opt) {
	  case 'f':
		input_freq = atoi(optarg);
		have_freq=true;
		break;
	  case 'i':
		have_invert=true;
	  case 'h':
	  default:
		usage();
		exit(0);
		break;
	}
	  
	  //printf("Header:\n%s", header);
	
    pf.filenum = 0;           // This is the crucial one to set to initialize things
    pf.rows_per_file = 200;  // Need to set this based on PSRFITS_MAXFILELEN

    // Now set values for our hdrinfo structure
    pf.hdr.scanlen = 86400; // in sec
    strcpy(pf.hdr.observer, "A. Eintein");
    strcpy(pf.hdr.telescope, "Effelsberg");
    strcpy(pf.hdr.obs_mode, "SEARCH");
    strcpy(pf.hdr.backend, "roach2");
    strcpy(pf.hdr.source, "J1745-2900");
    strcpy(pf.hdr.frontend, "C+_rcvr");
    strcpy(pf.hdr.project_id, "P84-16");

	char date_obs[64];
	if (ascii_header_get (header, "UTC_START", "%s", date_obs))
	  strcpy(pf.hdr.date_obs, date_obs);
    strcpy(pf.hdr.poln_type, "LIN");
    strcpy(pf.hdr.poln_order, "IQUV");
    strcpy(pf.hdr.track_mode, "TRACK");
    strcpy(pf.hdr.cal_mode, "OFF");
    strcpy(pf.hdr.feed_mode, "FA");

	double sampling_interval=0.0;
	if (ascii_header_get (header, "TSAMP", "%lf", &sampling_interval))
	  pf.hdr.dt = sampling_interval*1e-6; // IMPORTANT: TSAMP is the sampling period in microseconds
	
	double freq;
	if (ascii_header_get (header, "FREQ", "%lf", &freq))
	  pf.hdr.fctr = freq;
	if (have_freq) pf.hdr.fctr = input_freq;
	
	double bw;
	if (ascii_header_get (header, "BW", "%lf", &bw))
	  pf.hdr.BW = bw;

	int nchan;
	if (ascii_header_get (header, "NCHAN", "%d", &nchan))
	  pf.hdr.nchan = nchan;

	char MJD_start[64], *end;
	if (ascii_header_get (header, "MJD_START", "%s", MJD_start))
	  pf.hdr.MJD_epoch = strtold(MJD_start, &end);

	printf("Finished reading DADA header\n");
	
    pf.hdr.ra2000 = 266.41666667;
    dec2hms(pf.hdr.ra_str, pf.hdr.ra2000/15.0, 0);
    pf.hdr.dec2000 = -29.00833333;
    dec2hms(pf.hdr.dec_str, pf.hdr.dec2000, 1);

    pf.hdr.azimuth = 123.123;
    pf.hdr.zenith_ang = 23.0;
    pf.hdr.beam_FWHM = 0.25;
    pf.hdr.start_lst = 10000.0;
    pf.hdr.start_sec = 25000.82736876;
    pf.hdr.start_day = 55000;
    pf.hdr.scan_number = 3;
    pf.hdr.rcvr_polns = 2;
    pf.hdr.onlyI = 0;
    pf.hdr.summed_polns = 0;
    pf.hdr.offset_subint = 0;
    pf.hdr.orig_nchan = pf.hdr.nchan;
    pf.hdr.orig_df = pf.hdr.df = pf.hdr.BW / pf.hdr.nchan;
    pf.hdr.nbits = 8;
    pf.hdr.npol = 4;
    pf.hdr.chan_dm = 0.0;
    pf.hdr.fd_hand = 1;
    pf.hdr.fd_sang = 0;
    pf.hdr.fd_xyph = 0;    
    pf.hdr.be_phase = 1;
    pf.hdr.nsblk = 8192;
    pf.hdr.ds_time_fact = 1;
    pf.hdr.ds_freq_fact = 1;
    sprintf(pf.basefilename, "%s_%s_%05d_%4d", pf.hdr.backend, pf.hdr.source, (int) pf.hdr.MJD_epoch, (int) pf.hdr.fctr);

    psrfits_create(&pf);


    // Now set values for our subint structure
    pf.sub.tsubint = pf.hdr.nsblk * pf.hdr.dt;
    pf.tot_rows = 0.0;
    pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;
    pf.sub.lst = pf.hdr.start_lst;
    pf.sub.ra = pf.hdr.ra2000;
    pf.sub.dec = pf.hdr.dec2000;
    // GD slaEqgal(pf.hdr.ra2000*DEGTORAD, pf.hdr.dec2000*DEGTORAD, 
    //&pf.sub.glon, &pf.sub.glat);
    //pf.sub.glon *= RADTODEG;
    //pf.sub.glat *= RADTODEG;
    pf.sub.feed_ang = 0.0;
    pf.sub.pos_ang = 0.0;
    pf.sub.par_ang = 0.0;
    pf.sub.tel_az = pf.hdr.azimuth;
    pf.sub.tel_zen = pf.hdr.zenith_ang;
    pf.sub.bytes_per_subint = (pf.hdr.nbits * pf.hdr.nchan * 
                               pf.hdr.npol * pf.hdr.nsblk) / 8;
    pf.sub.FITS_typecode = TBYTE;  // 11 = byte

    // Create and initialize the subint arrays
    pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
    dtmp = pf.hdr.fctr - 0.5 * pf.hdr.BW + 0.5 * pf.hdr.df;
    for (ii = 0 ; ii < pf.hdr.nchan ; ii++) {
	  pf.sub.dat_freqs[ii] = dtmp + ii * pf.hdr.df;
	  pf.sub.dat_weights[ii] = 1.0;
    }
    pf.sub.dat_offsets = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
    pf.sub.dat_scales = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
    for (ii = 0 ; ii < pf.hdr.nchan * pf.hdr.npol ; ii++) {
	  pf.sub.dat_offsets[ii] = 0.0;
	  pf.sub.dat_scales[ii] = 1.0;
    }
 

    //fseek(pfi, DADA_HEADER_SIZE, SEEK_SET);

	uint8_t swap, *sptr;
	
    // This is what you would update for each time sample (likely just
    // adjusting the pointer to point to your data)
    pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);
	
    // Here is the real data-writing loop
    do {
      // Update the pf.sub entries here for each subint
      // as well as the pf.sub.data pointer
      
      fread(pf.sub.rawdata, pf.sub.bytes_per_subint, 1, pfi);
	  
	  // Flip the band if needed
	  // Only works for 8-bit values
	  if (have_invert) {
		for (jj=0; jj<pf.hdr.nsblk * pf.hdr.npol; jj++) {
		  sptr = (uint8_t *) &pf.sub.rawdata[jj*pf.hdr.nchan];
		  for (ii=0; ii<nchan/2; ii++) {
			swap = sptr[nchan-1-ii];
			sptr[nchan-1-ii] = sptr[ii];
			sptr[ii] = swap;
		  }
		}
	  }
		
	  pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;
	  
	  psrfits_write_subint(&pf);

	  if feof (pfi) 
	} while (pf.T < pf.hdr.scanlen && !feof (pfi) && !pf.status);
    
    // Close the last file and cleanup
    fits_close_file(pf.fptr, &(pf.status));
    free(pf.sub.dat_freqs);
    free(pf.sub.dat_weights);
    free(pf.sub.dat_offsets);
    free(pf.sub.dat_scales);
    free(pf.sub.data);

    printf("Done.  Wrote %d subints (%f sec) in %d files.  status = %d\n", 
           pf.tot_rows, pf.T, pf.filenum, pf.status);

    exit(0);
}
