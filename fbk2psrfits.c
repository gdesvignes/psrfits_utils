#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <stdbool.h>
#include "psrfits.h"
#include "fbk_header.h"

#ifndef DEGTORAD
#define DEGTORAD 0.017453292519943295769236907684886127134428718885417
#endif
#ifndef RADTODEG
#define RADTODEG 57.29577951308232087679815481410517033240547246656
#endif


void usage() {
  printf(
		             "Usage: fbk2psrfits [options] input_filename_base or filenames\n"
					 "Options:\n"
					 "  -h, --help               Print this\n"
					 "  -F, --frontend=FRONTEND  Frontend name, e.g. S36, S60, S45, S110\n"
					 "  -o, --off_chan           Correct the frequency by half a chan width\n"
					 "  -s, --srcname=NAME       Update source name\n "
					 "  -S HAve Stokes\n"
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
	{"frontend",   1, NULL, 'F'},
	{"offset",   0, NULL, 'o'},
	{"srcname",   1, NULL, 's'},
  };
  
    int ii;
    double dtmp;
    char filename[256], *adr;
    struct psrfits pf;
    int c1,c2;
    double c3;
	int opt, opti;
    int scan_id, beam_id;
	bool have_Stokes = false;
	bool have_chan_offset = false;
	bool have_frontend = false; char frontend[3];
	bool have_srcname = false; char srcname[24];

	while ((opt=getopt_long(argc,argv,"ob:t:j:i:F:s:p:P:F:CuU:SAqh",long_opts,&opti))!=-1) {
	  switch (opt) {
	  case 'F':
		have_frontend = true;
		strncpy(frontend, optarg, 3);
		break;
	  case 'o':
		have_chan_offset = true;
		break;
	  case 's':
		have_srcname = true;
		strncpy(srcname, optarg, 24);
		break;
	  case 'S':
		have_Stokes = true;
		break;
	  case 'h':
	  default:
		usage();
		exit(0);
		break;
	  }

	}
    
    // Only set the basefilename and not "filename"
    // Also, fptr will be set by psrfits_create_searchmode()

    FILE *pfi;

    strcpy(filename, argv[optind]);
    printf("Filename = %s\n", filename);

    // Read filterbank                                                                              
    if((pfi=fopen(filename,"r"))==NULL) {
      printf("Cannot open file %s\n", filename);
      exit(-2);
    } else {
      printf("Opening file %s\n", filename);
    }

    /* try to read the header */
    long long header_size, numsamps;
    header iheader;
    if (!(header_size=read_header(pfi,&iheader))) {
      printf("Error reading header\n");
      exit(-2);
    }
    //printf("header_size = %lld\n", header_size);

    // Get the number of samples                                                                              
    numsamps=nsamples(filename, header_size, &iheader);
    //printf("numsamps = %lld  ra=%lf  dec=%lf\n", numsamps, iheader.src_raj, iheader.src_dej);
    //fflush(stdout);
    
    pf.filenum = 0;           // This is the crucial one to set to initialize things
    pf.rows_per_file = 4000;  // Need to set this based on PSRFITS_MAXFILELEN

    // Now set values for our hdrinfo structure
    pf.hdr.scanlen = 3600; // in sec
    strcpy(pf.hdr.observer, "Observer");
    strcpy(pf.hdr.telescope, "Effelsberg");
    strcpy(pf.hdr.obs_mode, "SEARCH");
    strcpy(pf.hdr.backend, "PFFTS");
    strncpy(pf.hdr.source, iheader.source_name, 24);
    strcpy(pf.hdr.frontend, "LBand");
    strcpy(pf.hdr.project_id, "HTRUN");
    strcpy(pf.hdr.date_obs, "2010-01-01T05:15:30.000");
    strcpy(pf.hdr.poln_type, "LIN");
    strcpy(pf.hdr.poln_order, "IQUV");
    strcpy(pf.hdr.track_mode, "TRACK");
    strcpy(pf.hdr.cal_mode, "OFF");
    strcpy(pf.hdr.feed_mode, "FA");
    //printf("filename = %s\n", filename); fflush(stdout);
    if(strstr(filename,"/") != NULL) { adr=(char *)strrchr(filename,'/')+1; sprintf(pf.filename,"%s",adr); }
    else strncpy(pf.filename, filename, 80);
    //printf("filename = %s\n", pf.filename); fflush(stdout);
    sscanf(pf.filename, "%d_%d_%d_%dbit.fil", &scan_id, &c1, &beam_id, &c2);
      
    pf.hdr.dt = iheader.tsamp;
    pf.hdr.fctr = iheader.fch1+iheader.foff*iheader.nchans/2. - iheader.foff/2.;
	if (have_chan_offset)
	  pf.hdr.fctr += iheader.foff/2.;
	
    pf.hdr.BW = iheader.nchans * iheader.foff;
    angle_split(iheader.src_raj, &c1, &c2, &c3);
    sprintf(pf.hdr.ra_str, "%2.2d:%2.2d:%07.4lf", c1, c2, c3);
    
    angle_split(iheader.src_dej, &c1, &c2, &c3);
    sprintf(pf.hdr.dec_str, "%2.2d:%2.2d:%07.4lf", c1, c2, c3);
    pf.hdr.azimuth = iheader.az_start;
    pf.hdr.zenith_ang = iheader.za_start;
    pf.hdr.beam_FWHM = 0.25;
    pf.hdr.start_lst = 0.0;
    pf.hdr.start_sec = 0.0;
    pf.hdr.start_day = 55000;
    pf.hdr.scan_number = scan_id;
    pf.hdr.rcvr_polns = 2;
    pf.hdr.summed_polns = 0;
    pf.hdr.offset_subint = 0;
    pf.hdr.nchan = iheader.nchans;
    pf.hdr.orig_nchan = pf.hdr.nchan;
    pf.hdr.orig_df = pf.hdr.df = pf.hdr.BW / pf.hdr.nchan;
    pf.hdr.nbits = iheader.nbits;
    pf.hdr.npol = iheader.nifs;
	pf.hdr.onlyI = 0;
	if (pf.hdr.npol==1) pf.hdr.onlyI = 1;
    pf.hdr.chan_dm = 0.0;
    pf.hdr.fd_hand = 1;
    pf.hdr.fd_sang = 0;
    pf.hdr.fd_xyph = 0;
    pf.hdr.be_phase = 1;
    pf.hdr.ibeam = beam_id;
    pf.hdr.nsblk = 2048;

    pf.hdr.MJD_epoch = iheader.tstart;  // Note the "L" for long double
    pf.hdr.ds_time_fact = 1;
    pf.hdr.ds_freq_fact = 1;
    //sprintf(pf.basefilename, "%s_%04d_%02d", pf.hdr.backend, pf.hdr.scan_number, pf.hdr.ibeam);

    //psrfits_create(&pf);


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
	
    if (pf.hdr.nbits==32) pf.sub.FITS_typecode = TFLOAT;
	else pf.sub.FITS_typecode = TBYTE;  // 11 = byte

	//
	if (have_frontend) {
	  if (strncmp(frontend,"S36",3)==0) {
		strcpy(pf.hdr.backend, "PSRIX");
		strcpy(pf.hdr.frontend, "S36");
		strcpy(pf.hdr.poln_type, "CIRC");
		strcpy(pf.hdr.poln_order, "AABBCRCI");
		pf.hdr.fd_hand = -1; // To conform PA and Stokes V to PSR/IEEE 
	  }
	  else if (strncmp(frontend,"S60",3)==0) {
		strcpy(pf.hdr.backend, "PSRIX");
		strcpy(pf.hdr.frontend, "S60");
		strcpy(pf.hdr.poln_type, "CIRC");
		strcpy(pf.hdr.poln_order, "AABBCRCI");		
	  }
	  else if (strncmp(frontend,"S110",3)==0) {
		strcpy(pf.hdr.backend, "PSRIX");
		strcpy(pf.hdr.frontend, "S110");
		strcpy(pf.hdr.poln_type, "CIRC");
		strcpy(pf.hdr.poln_order, "AABBCRCI");
		pf.hdr.be_phase = -1; // To conform Stokes V to PSR/IEEE
	  }
	  else if (strncmp(frontend,"S45",3)==0) {
		strcpy(pf.hdr.backend, "PSRIX");
		strcpy(pf.hdr.frontend, "S45");
		strcpy(pf.hdr.poln_type, "LIN");
		strcpy(pf.hdr.poln_order, "AABBCRCI");
		pf.hdr.be_phase = -1; // To conform Stokes V to PSR/IEEE
	  }
	  else {printf("Frontend %s not recognised\n", frontend);}
	}
	if (have_srcname) strncpy(pf.hdr.source, srcname, 24);

	if (have_Stokes) strcpy(pf.hdr.poln_order, "IQUV");
	
	if (strncmp(pf.hdr.backend, "PFFTS",5)==0)
	  sprintf(pf.basefilename, "%s_%04d_%02d", pf.hdr.backend, pf.hdr.scan_number, pf.hdr.ibeam);
	else
	  sprintf(pf.basefilename, "%s_%s_%5d_%d", pf.hdr.backend, pf.hdr.source, (int) pf.hdr.MJD_epoch, (int) pf.hdr.fctr);
	psrfits_create(&pf);
	
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
 

    // This is what you would update for each time sample (likely just
    // adjusting the pointer to point to your data)
    pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);
   
    // Here is the real data-writing loop
    do {
      // Update the pf.sub entries here for each subint
      // as well as the pf.sub.data pointer
      
      fread(pf.sub.rawdata, pf.sub.bytes_per_subint, 1, pfi);
      pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;
      
      psrfits_write_subint(&pf);
    } while (pf.T < pf.hdr.scanlen && !feof (pfi) && !pf.status);
    
    // Close the last file and cleanup
    fits_close_file(pf.fptr, &(pf.status));
    free(pf.sub.dat_freqs);
    free(pf.sub.dat_weights);
    free(pf.sub.dat_offsets);
    free(pf.sub.dat_scales);
    free(pf.sub.rawdata);

    printf("Done.  Wrote %d subints (%f sec) in %d files.  status = %d\n", 
           pf.tot_rows, pf.T, pf.filenum, pf.status);

    exit(0);
}
