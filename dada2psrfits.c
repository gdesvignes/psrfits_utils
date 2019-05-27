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
		 "  -c, --flipV              Flip polarisation 4, e.g. Stokes V\n"
		 "  -h, --help               Print this\n"
		 "  -f fff, --freq=ffff      Set the centre frequency (MHz)\n"
		 "  -i, --invert             Invert band\n"
		 "  -p, --psrname            Set PSRNAME string in psrfits archive\n"
		 "  -r, --ra                 Set RA string in psrfits archive\n"
		 "  -d, --dec                Set DEC string in psrfits archive\n"
		 "  -v, --pv                 Set Pico Veleta as the observatory\n"
		 "  -s, --Slin2circ          Transform from linear to circular feed Stokes output\n"
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
	{"flipV",   0, NULL, 'c'},
	{"freq",    1, NULL, 'f'},
	{"invert",  0, NULL, 'i'},
	{"help",    0, NULL, 'h'},
	{"psrname", 0, NULL, 'p'},
	{"ra",      0, NULL, 'r'},
	{"dec",     0, NULL, 'd'},
	{"pv",      0, NULL, 'v'},
	{"Slin2circ",0, NULL, 's'},
	{0,0,0,0}
  };

  int opt, opti;
    int ii, jj;
    double dtmp;
    struct psrfits pf;
    long long nspec;    
	double input_freq;
	bool have_freq=false, have_invert=false, do_break=false, is_pico=false;
	bool have_ra=false, have_dec=false, have_source=false, have_flipV=false;
	bool have_lin2circ=false;
	char ra_str[16], dec_str[16], source[24];
	while ((opt=getopt_long(argc,argv,"cd:f:ir:p:svh", long_opts,&opti))!=-1) {
	  switch (opt) {
	  case 'c':
		have_flipV=true;
		break;
	  case 'd':
		strncpy(dec_str, optarg, 16);
		have_dec=true;
		break;
	  case 'f':
		input_freq = atoi(optarg);
		have_freq=true;
		break;
	  case 'i':
		have_invert=true;
		break;
	  case 'r':
		strncpy(ra_str, optarg, 16);
		have_ra=true;
		break;
	  case 'p':
		strncpy(source, optarg, 24);
		have_source=true;
		break;
	  case 'v':
		is_pico=true;
		break;
	  case 's':
		have_lin2circ=true;
	  case 'h':
	  default:
		usage();
		exit(0);
		break;
	  }
	}

	if (optind==argc) {
	  usage();
	  exit(1);
	}

    FILE *pfi;
	char filename[256], sfn[256];
	strcpy(filename, argv[optind]);
	uint64_t fs;
	//sscanf(filename,"%*[^_]_%"PRIu64".000000.dada", &fs);
	sscanf(filename,"%[^_]_%"PRIu64".000000.dada", sfn, &fs);
	printf("sfn = %s  fs = %"PRIu64"\n", sfn, fs);
	printf("Opening file %s\n", filename);
	pfi = fopen(filename, "r");

	char header[DADA_HEADER_SIZE];

	if (fread (header, 1, DADA_HEADER_SIZE, pfi) != DADA_HEADER_SIZE) {
	  fprintf (stderr, "Error getting header from %s\n", filename);
	  exit(-1);
	}
	
    pf.filenum = 0;           // This is the crucial one to set to initialize things
    pf.rows_per_file = 200;  // Need to set this based on PSR
	pf.hdr.nbits = 8;
    // Now set values for our hdrinfo structure
    pf.hdr.scanlen = 86400; // in sec
    strcpy(pf.hdr.observer, "A. Eintein");
	if (is_pico) {
	  strcpy(pf.hdr.telescope, "Pico Veleta");
	  pf.hdr.be_phase = 1;
	} else {
	  strcpy(pf.hdr.telescope, "Effelsberg");
	  pf.hdr.be_phase = -1;
	}
    strcpy(pf.hdr.obs_mode, "SEARCH");
    strcpy(pf.hdr.backend, "roach2");
    strcpy(pf.hdr.source, source);
    strcpy(pf.hdr.frontend, "S45-1");
    strcpy(pf.hdr.project_id, "GC");

	char date_obs[64];
	if (ascii_header_get (header, "UTC_START", "%s", date_obs))
	  strcpy(pf.hdr.date_obs, date_obs);
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
	if (have_invert) pf.hdr.BW *= -1;

	int nchan;
	if (ascii_header_get (header, "NCHAN", "%d", &nchan))
	  pf.hdr.nchan = nchan;

	int nbits;
	if (ascii_header_get (header, "NBIT", "%d", &nbits))
	  pf.hdr.nbits = nbits;

	char basis[16];
	if (ascii_header_get (header, "BASIS", "%s", basis)) {
	  if (strncmp("Circular", basis, 8)==0) { strcpy(pf.hdr.poln_type, "CIRC");}
	  else { strcpy(pf.hdr.poln_type, "LIN");}
	} else strcpy(pf.hdr.poln_type, "LIN");

	//char source[24];
	if (have_source)
	  strncpy(pf.hdr.source, source, 24);
	else if (ascii_header_get (header, "SOURCE", "%s", &source))
	  strncpy(pf.hdr.source, source, 24);

	//char ra_str[16];
	if (have_ra)
	  strncpy(pf.hdr.ra_str, ra_str, 16);
	else if (ascii_header_get (header, "RA", "%s", &ra_str))
	  strncpy(pf.hdr.ra_str, ra_str, 16);
	
	//char dec_str[16];
	if (have_dec)
	  strncpy(pf.hdr.dec_str, dec_str, 16);
	else if (ascii_header_get (header, "DEC", "%s", &dec_str))
	  strncpy(pf.hdr.dec_str, dec_str, 16);

	char MJD_start[64], *end;
	if (ascii_header_get (header, "MJD_START", "%s", MJD_start))
	  pf.hdr.MJD_epoch = strtold(MJD_start, &end);

	uint64_t file_size;
	ascii_header_get (header, "FILE_SIZE", "%"PRIu64"", &file_size);

	uint64_t bytes_per_second;
	ascii_header_get (header, "BYTES_PER_SECOND", "%"PRIu64"", &bytes_per_second);
	
	uint64_t offset;
	ascii_header_get (header, "OBS_OFFSET", "%"PRIu64"", &offset);
	if (offset) {
	  long double time_offset = 0.;
	  time_offset = offset / bytes_per_second;
	  printf("Initial file has OBS_OFFSET!=0. Adding %lf second to MJD_epoch\n", time_offset);
	}

	// TODO: should match with telescope Effelsberg
	if (pf.hdr.MJD_epoch < 58472 && pf.hdr.fctr == 7000.) {
	  printf("Correcting for the wrong BW sign in early EffelsbergPSRIX2 data\n");
	  pf.hdr.BW *= -1;
	}
	
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
    pf.hdr.npol = 4;
    pf.hdr.chan_dm = 0.0;
    pf.hdr.fd_hand = 1;
    pf.hdr.fd_sang = 0;
    pf.hdr.fd_xyph = 0;    
    pf.hdr.nsblk = 8192;
    pf.hdr.ds_time_fact = 1;
    pf.hdr.ds_freq_fact = 1;

	sprintf(pf.hdr.ra_str, ra_str);
	sprintf(pf.hdr.dec_str, dec_str);

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
	if (pf.hdr.nbits==32) {pf.sub.FITS_typecode = TFLOAT;
	  printf("Input is 32b. Output FITS type defaulting to FLOAT\n");}
    else pf.sub.FITS_typecode = TBYTE;  // 11 = byte

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
 

	uint8_t swap, *uptr;
	int8_t *sptr;
	char *buf_L2C;
	int ipol, opol;
    // This is what you would update for each time sample (likely just
    // adjusting the pointer to point to your data)
    pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);
	buf_L2C = (unsigned char *)malloc(pf.hdr.nchan*pf.hdr.nbits/8);
	
    // Here is the real data-writing loop
    do {
      // Update the pf.sub entries here for each subint
      // as well as the pf.sub.data pointer
      nspec = fread(pf.sub.rawdata, pf.hdr.nbits/8 * pf.hdr.nchan * pf.hdr.npol, pf.hdr.nsblk, pfi);
	  
      if (feof (pfi) ) {
		// Close file
		printf("Closing file %s\n", filename);
		fclose(pfi);
		fs += file_size;
		sprintf(filename, "%s_%016"PRIu64".000000.dada", sfn, fs);

		if ( (pfi = fopen(filename, "r")) != NULL) {
		  printf("Opening file %s\n", filename);
		  fread (header, 1, DADA_HEADER_SIZE, pfi);
		}
		else {
		  break;
		}
	  }

	  // Need to read some more data from new file if there is 
	  if (nspec < pf.hdr.nsblk) {
	    fread(&pf.sub.rawdata[nspec*pf.hdr.nbits/8 * pf.hdr.nchan * pf.hdr.npol],
					  pf.hdr.nbits * pf.hdr.nchan * pf.hdr.npol/8, pf.hdr.nsblk-nspec, pfi);
		//nspec+=nspec2;

	  }

	  if (have_lin2circ && (strncmp(pf.hdr.poln_type, "CIRC",4)==0)) {
		for (jj=0; jj<pf.hdr.nsblk;jj++) {
		  // First copy Stokes Q to tmp
		  ipol = 1;
		  memcpy(buf_L2C, &pf.sub.rawdata[jj*pf.hdr.npol * pf.hdr.nchan * ipol * pf.hdr.nchan], pf.hdr.nchan * pf.hdr.nbits/8);
		  // Copy Stokes U to Q
		  opol = 1;
		  ipol = 2;
		  memcpy(&pf.sub.rawdata[jj*pf.hdr.npol * pf.hdr.nchan + opol * pf.hdr.nchan],
				 &pf.sub.rawdata[jj*pf.hdr.npol * pf.hdr.nchan + ipol * pf.hdr.nchan], pf.hdr.nchan * pf.hdr.nbits/8);
		  // Copy Stokes V to U
		  opol = 2;
		  ipol = 3;
		  memcpy(&pf.sub.rawdata[jj*pf.hdr.npol * pf.hdr.nchan + opol * pf.hdr.nchan],
				 &pf.sub.rawdata[jj*pf.hdr.npol * pf.hdr.nchan + ipol * pf.hdr.nchan], pf.hdr.nchan * pf.hdr.nbits/8);
		  // Copy tmp to Stokes V
		  opol = 3;
		  memcpy(&pf.sub.rawdata[jj*pf.hdr.npol * pf.hdr.nchan + opol * pf.hdr.nchan], buf_L2C, pf.hdr.nchan * pf.hdr.nbits/8);
		}
	  }
	  
	  
	  // Flip the band if needed
	  // Only works for 8-bit values
	  if (have_invert) {
		if (pf.hdr.nbits!=8) {printf("!8bits band inversion not supported yet. Exit\n" );return(-1);}
		for (jj=0; jj<pf.hdr.nsblk * pf.hdr.npol; jj++) {
		  uptr = (uint8_t *) &pf.sub.rawdata[jj*pf.hdr.nchan];
		  for (ii=0; ii<nchan/2; ii++) {
			swap = uptr[nchan-1-ii];
			uptr[nchan-1-ii] = uptr[ii];
			uptr[ii] = swap;
		  }
		}
	  }

	  if (have_flipV) {
		if (pf.hdr.nbits!=8) {printf("!8bits flip Stokes V not supported yet. Exit\n" );return(-1);}
		for (jj=0; jj<pf.hdr.nsblk; jj++) {
		  sptr = (int8_t *) &pf.sub.rawdata[jj*pf.hdr.nchan*pf.hdr.npol + 3*pf.hdr.nchan];
		  for (ii=0; ii<pf.hdr.nchan; ii++) {
			*sptr *= -1;
			sptr++;
		  }
		}
	  }
		
	  pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;
	  
	  psrfits_write_subint(&pf);

	  //} while (pf.T < pf.hdr.scanlen && !feof (pfi) && !pf.status);
	} while (pf.T < pf.hdr.scanlen && !feof (pfi));
    
    // Close the last file and cleanup
    fits_close_file(pf.fptr, &(pf.status));
    free(pf.sub.dat_freqs);
    free(pf.sub.dat_weights);
    free(pf.sub.dat_offsets);
    free(pf.sub.dat_scales);
    free(pf.sub.data);
	free(buf_L2C);
    printf("Done.  Wrote %d subints (%f sec) in %d files.  status = %d\n", 
           pf.tot_rows, pf.T, pf.filenum, pf.status);

    exit(0);
}
