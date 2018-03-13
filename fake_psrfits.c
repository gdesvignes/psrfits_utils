#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <getopt.h>
#include <stdbool.h>
#include <inttypes.h>
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
		 "Usage: dada2psrfits [options] input_filename_base\n"
		 "Options:\n"
		 "  -b, --b2020              Set default options for B2020+28 observations\n"
		 "  -h, --help               Print this\n"
		 "  -f fff, --freq=ffff      Set the centre frequency (MHz)\n"
		 "  -i, --invert             Invert band\n"
		 "  -j, --j2205              Set default options for J2205+6012 observations\n"
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
	{"b2020",   0, NULL, 'b'},
	{"j2205",   0, NULL, 'j'},
	{"freq",    1, NULL, 'f'},
	{"invert",  0, NULL, 'i'},
	{"help",    0, NULL, 'h'},
	{0,0,0,0}
  };

  int opt, opti;
  int ii, jj, status;
    double dtmp;
    struct psrfits pf;
    long long nspec;    
	double input_freq;
	bool have_freq=false, have_invert=false, is_b2020=false, do_break=false, is_j2205=false;
	while ((opt=getopt_long(argc,argv,"bjf:ih", long_opts,&opti))!=-1) {
	  switch (opt) {
	  case 'b':
		is_b2020=true;
		break;
	  case 'j':
		is_j2205=true;
		break;
	  case 'f':
		input_freq = atoi(optarg);
		have_freq=true;
		break;
	  case 'i':
		have_invert=true;
		break;
	  case 'h':
	  default:
		usage();
		exit(0);
		break;
	  }
	}


    FILE *pfi;
	char filename[256], sfn[256];
	uint64_t fs;
	//sscanf(filename,"%*[^_]_%"PRIu64".000000.dada", &fs);

	char header[DADA_HEADER_SIZE];
	
    pf.filenum = 0;           // This is the crucial one to set to initialize things
    pf.rows_per_file = 200;  // Need to set this based on PSR

    // Now set values for our hdrinfo structure
    pf.hdr.scanlen = 600; // in sec
    strcpy(pf.hdr.observer, "A. Eintein");
    strcpy(pf.hdr.telescope, "Effelsberg");
    strcpy(pf.hdr.obs_mode, "SEARCH");
    strcpy(pf.hdr.backend, "roach2");
    strcpy(pf.hdr.source, "J1745-2900");
    strcpy(pf.hdr.frontend, "C+_rcvr");
    strcpy(pf.hdr.project_id, "P84-16");

	char date_obs[64];

	strcpy(pf.hdr.date_obs, date_obs);
    strcpy(pf.hdr.poln_type, "LIN");
    strcpy(pf.hdr.poln_order, "IQUV");
    strcpy(pf.hdr.track_mode, "TRACK");
    strcpy(pf.hdr.cal_mode, "OFF");
    strcpy(pf.hdr.feed_mode, "FA");

	double sampling_interval=64.0;
	pf.hdr.dt = sampling_interval*1e-6; // IMPORTANT: TSAMP is the sampling period in microseconds
	
	double freq=1512;
	pf.hdr.fctr = freq;
	
	double bw=512;
	pf.hdr.BW = bw;

	int nchan=512;
	pf.hdr.nchan = nchan;
	  
	char MJD_start[64], *end;

	uint64_t file_size;

	pf.hdr.MJD_epoch = 55555.123123123123123123L;
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

	if (is_b2020) {
	  sprintf(pf.hdr.source, "B2020+28");
	  sprintf(pf.hdr.ra_str, "20:22:37.070");
	  sprintf(pf.hdr.dec_str, "+28:54:23.10");
	}

	if (is_j2205) {
	  sprintf(pf.hdr.source, "J2205+6012");
	  sprintf(pf.hdr.ra_str, "22:05:34.2016");
	  sprintf(pf.hdr.dec_str, "+60:12:55.1417");
	}
	
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
 

	uint8_t swap, *sptr;
	
    // This is what you would update for each time sample (likely just
    // adjusting the pointer to point to your data)
    pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);
	int percolnum;
	
    // Here is the real data-writing loop
    do {
      // Update the pf.sub entries here for each subint
      // as well as the pf.sub.data pointer
      //nspec = fread(pf.sub.rawdata, pf.hdr.nbits * pf.hdr.nchan * pf.hdr.npol/8, pf.hdr.nsblk, pfi);


	  // Need to read some more data from new file if there is 

		
	  pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;
	  //printf("offs = %lf\n", pf.sub.offs);
	  psrfits_write_subint(&pf);
	  fits_get_colnum(pf.fptr, CASEINSEN, "PERIOD", &percolnum, &status);
	  fits_get_colnum(pf.fptr, CASEINSEN, "OFFS_SUB", &percolnum, &status);
	  fits_get_colnum(pf.fptr, CASEINSEN, "DAT_FREQ", &percolnum, &status);
	  fits_get_colnum(pf.fptr, CASEINSEN, "TSUBINT", &percolnum, &status);
	  
	  //} while (pf.T < pf.hdr.scanlen && !feof (pfi) && !pf.status);
	} while (pf.T < pf.hdr.scanlen);
    
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
