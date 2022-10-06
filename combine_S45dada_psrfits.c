// Only works so far for 8-b and 32-b DADA in Stokes Search mode
// G. Desvignes, April 2019

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

#define ZAPLIST "/fpra/galc/01/gdesvignes/src/psrfits_utils/S45_channels.zaplist"

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



void usage() {
  printf(
		 "Usage: dada2psrfits [options] dadafile1 dadafile2\n"
		 "Options:\n"
		 "  -h, --help               Print this\n"
		 "  -f fff, --freq=ffff      Set the centre frequency (MHz)\n"
		 "  -p, --psrname            Set PSRNAME string in psrfits archive\n"
		 "  -r, --ra                 Set RA string in psrfits archive (e.g. 18:09:51.090) \n"
		 "  -d, --dec                Set DEC string in psrfits archive (e.g. -19:43:51.90)\n"
		 "  -v, --pv                 Set Pico Veleta as the observatory\n"
		 "  -z, --zap                Use zap list for Effelsberg S45 data\n"
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
	{"help",    0, NULL, 'h'},
	{"psrname", 0, NULL, 'p'},
	{"ra",      0, NULL, 'r'},
	{"dec",     0, NULL, 'd'},
	{"pv",      0, NULL, 'v'},
	{"zap",     0, NULL, 'z'},
	{0,0,0,0}
  };

  int opt, opti;
  int ii, jj;
  double dtmp;
  struct psrfits pf;
  long long nspec, nspec1, nspec2;    
  double input_freq;
  bool have_freq=false, do_break=false, is_pico=false;
  bool have_ra=false, have_dec=false, have_source=false, have_zap=false;
  char ra_str[16], dec_str[16], source[24];
  while ((opt=getopt_long(argc,argv,"d:f:r:p:vzh", long_opts,&opti))!=-1) {
      switch (opt) {
       case 'd':
	  strncpy(dec_str, optarg, 16);
	  have_dec=true;
	  break;
      case 'f':
	  input_freq = atoi(optarg);
	  have_freq=true;
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
	  case 'z':
		have_zap=true;
		break;
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
  
  FILE *pfi, *pfi2;
  char *adr, filename[256], sfn[256], filename2[256], sfn2[256], path[256], path2[256];
  char ffn[256], ffn2[256];
  char sfx[40];
  strcpy(ffn, argv[optind]);
  if((adr=(char *)strrchr(ffn,'/')) != NULL) strcpy(filename,adr+1);
  strcpy(ffn2, argv[optind+1]);
  if((adr=(char *)strrchr(ffn2,'/')) != NULL) strcpy(filename2,adr+1);
  
  strcpy(path, ffn);
  strcpy(path2, ffn2);
  if((adr=(char *)strrchr(path,'/')) != NULL) sprintf(adr,"\0");
  if((adr=(char *)strrchr(path2,'/')) != NULL) sprintf(adr,"\0");
  

  uint64_t fs, fs2;
  sscanf(filename,"%[^_]_%"PRIu64".%s", sfn, &fs, sfx);
  printf("Opening %s\n", ffn);
  sscanf(filename2,"%[^_]_%"PRIu64".%s", sfn2, &fs2, sfx);
  printf("Opening %s\n", ffn2);

  if (strncmp(sfx, "R2D", 3)==0) {
      printf("New R2D format\n");
  } else {
      printf("Old format: %s\n", sfx);
  }
  
  if (fs!=fs2) {printf("Error with input files:\n%s != %s\n", ffn, ffn2); return(-1);}

  pfi = fopen(ffn, "r");
  pfi2 = fopen(ffn2, "r");


  if (pfi==NULL) {
	printf("Can not open %s\n", ffn);
	return -1;
  }
  if (pfi2==NULL) {
	printf("Can not open %s\n", ffn2);
	return -1;
  }
  
  char header[DADA_HEADER_SIZE], header2[DADA_HEADER_SIZE];  
  if (fread (header, 1, DADA_HEADER_SIZE, pfi) != DADA_HEADER_SIZE) {
      fprintf (stderr, "Error getting header from %s\n", ffn);
      exit(-1);
  }

  if (fread (header2, 1, DADA_HEADER_SIZE, pfi2) != DADA_HEADER_SIZE) {
      fprintf (stderr, "Error getting header from %s\n", ffn2);
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
	  strcpy(pf.hdr.frontend, "PV-RCVR");
      pf.hdr.be_phase = 1;
  } else {
      strcpy(pf.hdr.telescope, "Effelsberg");
	  strcpy(pf.hdr.frontend, "S45-1");
      pf.hdr.be_phase = -1;
  }
  strcpy(pf.hdr.obs_mode, "SEARCH");
  strcpy(pf.hdr.backend, "roach2");
  strcpy(pf.hdr.source, source);
  strcpy(pf.hdr.project_id, "GC");

  char date_obs[64], date_obs2[64];
  if (ascii_header_get (header, "UTC_START", "%s", date_obs) && ascii_header_get (header2, "UTC_START", "%s", date_obs2)) {
      strcpy(pf.hdr.date_obs, date_obs);  
  } else {
      printf("Error: UTC_START Band 1: %s != Band 2: %s\n", date_obs, date_obs2);
	  exit(-1);
  } 

  strcpy(pf.hdr.poln_type, "LIN");
  strcpy(pf.hdr.poln_order, "IQUV");
  strcpy(pf.hdr.track_mode, "TRACK");
  strcpy(pf.hdr.cal_mode, "OFF");
  strcpy(pf.hdr.feed_mode, "FA");

  double sampling_interval=0.0;
  if (ascii_header_get (header, "TSAMP", "%lf", &sampling_interval))
      pf.hdr.dt = sampling_interval*1e-6; // IMPORTANT: TSAMP is the sampling period in microseconds
	
  double freq, freq2;
  if (ascii_header_get (header, "FREQ", "%lf", &freq) && ascii_header_get (header2, "FREQ", "%lf", &freq2)) {
	if (freq>freq2) {
	  printf("Error: Freq file 1: %lf  > Freq file 2: %lf\nSwap filenames in command line!\n", freq, freq2);
	  //exit(-1);
	}
      pf.hdr.fctr = (freq + freq2) / 2.;
  }
  if (have_freq) pf.hdr.fctr = input_freq;
  printf("New central frequency: %lf\n", pf.hdr.fctr);
  
  double bw, bw2;
  if (ascii_header_get (header, "BW", "%lf", &bw) && ascii_header_get (header2, "BW", "%lf", &bw2))
      pf.hdr.BW = 2*bw; // Assumes sign of BW from band 1
  //if (have_invert) pf.hdr.BW *= -1;

  int nchan;
  if (ascii_header_get (header, "NCHAN", "%d", &nchan))
      pf.hdr.nchan = 2*nchan;
  
  int nbits;
  if (ascii_header_get (header, "NBIT", "%d", &nbits))
      pf.hdr.nbits = nbits;
  
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
  
  char MJD_start[64],  MJD_start2[64], *end;
  if (ascii_header_get (header, "MJD_START", "%s", MJD_start) && ascii_header_get (header2, "MJD_START", "%s", MJD_start2)) {
      if (strcmp(MJD_start, MJD_start2))
	  printf("Error: MJD_START band 1 != band2\n");
      else pf.hdr.MJD_epoch = strtold(MJD_start, &end);
  }
  
  uint64_t file_size;
  ascii_header_get (header, "FILE_SIZE", "%"PRIu64"", &file_size);
  file_size = fs; // GD: File_size can no longer be trusted. Recompute file_size manually while reading the data
  // If we provide the first file, fs should be zero.
  
  uint64_t bytes_per_second;
  ascii_header_get (header, "BYTES_PER_SECOND", "%"PRIu64"", &bytes_per_second);
  
  uint64_t offset;
  ascii_header_get (header, "OBS_OFFSET", "%"PRIu64"", &offset);
  if (offset) {
      long double time_offset = 0.;
      time_offset = (long double) offset / (long double) bytes_per_second;
      printf("Initial file has OBS_OFFSET=%llf, bytes_per_second=%llf. Adding %llf second to MJD_epoch\n", (long double) offset, (long double) bytes_per_second, time_offset);
      pf.hdr.MJD_epoch += time_offset/86400;
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

  if (have_zap) {
	// Read zap list
	FILE *pfzap;
	char line[8];
	int nchan2zap=0;
	int chan2zap[2048];
	if ((pfzap=fopen(ZAPLIST, "r"))==NULL) {
	  printf("Could not open zaplist %s\n", ZAPLIST);
	  exit(-1);
	} else {
	  printf("Opening zaplist %s\n", ZAPLIST);
	  while(!feof(pfzap) && fgets(line,8,pfzap)) {
		sscanf(line, "%d", &chan2zap[nchan2zap]);
		nchan2zap++;
	  }
	  fclose(pfzap);
	}
	printf("%d channels to zap\n", nchan2zap);

	// Apply the weights
	for (ii=0; ii<nchan2zap; ii++) {
            if (chan2zap[ii]<0 || chan2zap[ii] > pf.hdr.nchan) {
                printf("Channel %d outside of channels range in PSRFITS file 0-%d. Skipping.\n", chan2zap[ii], pf.hdr.nchan);
                continue;
            }
            pf.sub.dat_weights[chan2zap[ii]] = 0.0;
	}
  }      
	

  uint8_t swap, *uptr;
  int8_t *sptr;
  char *tmpdata, *tmpdata2;
  float *tmpdata2_f, *fptr;
  tmpdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);
  tmpdata2 = &tmpdata[pf.sub.bytes_per_subint/2]; 
  tmpdata2_f = (float *) &tmpdata[pf.sub.bytes_per_subint/2];
  pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);

  // Here is the real data-writing loop
  do {
      // Read one subint from each of the two files
      nspec1 = fread(tmpdata, pf.hdr.nbits/8 * nchan * pf.hdr.npol, pf.hdr.nsblk, pfi);
      nspec2 = fread(tmpdata2, pf.hdr.nbits/8 * nchan * pf.hdr.npol, pf.hdr.nsblk, pfi2);

      // Doesn't matter, EOF check on both files would stop the program
      nspec = nspec1;
      file_size += nspec * pf.hdr.nbits/8 * nchan * pf.hdr.npol; // Assume we read the same number of spectra
      
      // Check if EOF and open next file from band1
      if (feof (pfi)) {
	  // Close file
	  printf("Closing file %s\n", filename);
	  fclose(pfi);
	  sprintf(filename, "%s/%s_%016"PRIu64".%s", path, sfn, file_size, sfx);
	  
	  if ( (pfi = fopen(filename, "r")) != NULL) {
	      printf("Opening file %s\n", filename);
	      fread (header, 1, DADA_HEADER_SIZE, pfi);
	  }
	  else {
	      printf("Trying to open file %s\n", filename);
	      printf("No more data. Finishing...\n");
	      break;
	  }
      }

      // Check if EOF and open next file from band2
      if (feof (pfi2) ) {
	  // Close file
	  printf("Closing file %s\n", filename2);
          fclose(pfi2);
          sprintf(filename2, "%s/%s_%016"PRIu64".%s", path2, sfn2, file_size, sfx);

          if ( (pfi2 = fopen(filename2, "r")) != NULL) {
              printf("Opening file %s\n", filename2);
              fread (header2, 1, DADA_HEADER_SIZE, pfi2);
          }
          else {
              break;
          }
      }

      // Need to read some more data from both files
      if (nspec < pf.hdr.nsblk) {
		fread(&tmpdata[nspec*pf.hdr.nbits/8 * nchan * pf.hdr.npol],
			  pf.hdr.nbits * nchan * pf.hdr.npol/8, pf.hdr.nsblk-nspec, pfi);
		fread(&tmpdata2[nspec*pf.hdr.nbits/8 * nchan * pf.hdr.npol],
			  pf.hdr.nbits * nchan * pf.hdr.npol/8, pf.hdr.nsblk-nspec, pfi2);
		file_size += (pf.hdr.nsblk-nspec) * pf.hdr.nbits/8 * nchan * pf.hdr.npol;
      }
	  
      // Flip the second band
      // Only works for 8-bit values
      if (pf.hdr.nbits==8) {
		for (jj=0; jj<pf.hdr.nsblk; jj++) {
		  for (ii=0; ii< 3; ii++) {
			uptr = (uint8_t *) &tmpdata2[(jj*pf.hdr.npol + ii) * nchan];
			for (int kk=0; kk<nchan/2; kk++) {
			  swap = uptr[nchan-1-kk];
			  uptr[nchan-1-kk] = uptr[kk];
			  uptr[kk] = swap;
			}
		  }
		  // Flip the sign of 4th polarisation, Stokes V, while flipping the band
		  for (ii=3; ii<4; ii++) {
			sptr = (int8_t *) &tmpdata2[(jj*pf.hdr.npol + ii) * nchan];
			for (int kk=0; kk<nchan/2; kk++) {
			  swap = -1*sptr[nchan-1-kk];
			  sptr[nchan-1-kk] = -1*sptr[kk];
			  sptr[kk] = swap;
			}
		  }
		}
	  }

	  // Flip the second band
	  // Only works for 32-bit values
	  if (pf.hdr.nbits==32) {
		for (jj=0; jj<pf.hdr.nsblk; jj++) {
		  for (ii=0; ii< 3; ii++) {
			fptr = (float *) &tmpdata2_f[(jj*pf.hdr.npol + ii) * nchan];
			for (int kk=0; kk<nchan/2; kk++) {
			  swap = fptr[nchan-1-kk];
			  fptr[nchan-1-kk] = fptr[kk];
			  fptr[kk] = swap;
			}
		  }
		  for (ii=3; ii<4; ii++) {
			fptr = (float *) &tmpdata2_f[(jj*pf.hdr.npol + ii) * nchan];
			for (int kk=0; kk<nchan/2; kk++) {
			  swap = -1*fptr[nchan-1-kk];
			  fptr[nchan-1-kk] = -1*fptr[kk];
			  fptr[kk] = swap;
			}
		  }
		}
	  }

      // Reorder the data / Transpose from the two separated bands to one band
      reorder_data(pf.sub.rawdata, tmpdata, 2, pf.hdr.nsblk, pf.hdr.npol, nchan, pf.hdr.nbits);
      
      pf.sub.offs = (pf.tot_rows + 0.5) * pf.sub.tsubint;
      
      psrfits_write_subint(&pf);
      
      //} while (pf.T < pf.hdr.scanlen && !feof (pfi) && !pf.status);
  } while (pf.T < pf.hdr.scanlen && !feof (pfi)  && !feof (pfi2));
  
  // Close the last file and cleanup
  fits_close_file(pf.fptr, &(pf.status));
  free(pf.sub.dat_freqs);
  free(pf.sub.dat_weights);
  free(pf.sub.dat_offsets);
  free(pf.sub.dat_scales);
  free(pf.sub.rawdata);
  free(tmpdata);
  printf("Done.  Wrote %d subints (%f sec) in %d files.  status = %d\n", 
	 pf.tot_rows, pf.T, pf.filenum, pf.status);
  
  exit(0);
}
