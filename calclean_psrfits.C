#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <getopt.h>
#include <stdbool.h>
#include <inttypes.h>
#include "psrfits.h"

#ifdef HAVE_PSRCHIVE
#include <PolnCalibrator.h>
#endif

void usage() {
  printf(
		 "Usage: dada2psrfits [options] input_filename_base\n"
		 "Options:\n"
		 "  -h, --help               Print this\n"
		 );
}

int main(int argc, char *argv[]) {

  /* Cmd line */
  static struct option long_opts[] = {
	{"help",    0, NULL, 'h'},
	{"cal",     1, NULL, 'c'},
	{0,0,0,0}
  };

  int status,opt, opti;
  struct psrfits pfi, pf;
  char ra_str[16], dec_str[16], source[24];
  char cal_file[128];
  bool have_cal_file=false;
  while ((opt=getopt_long(argc,argv,"c:h", long_opts,&opti))!=-1) {
      switch (opt) {
      case 'c':
	  strncpy(cal_file, optarg, 128);
	  have_cal_file = true;
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

  // Open input file
  status = psrfits_open(&pfi);
  if (status) fits_report_error(stderr, status);

  // Change to 32 bits
  pf.hdr.nbits = 32;

  // Alloc data buffers for the output PSRFITS files
  pf.sub.bytes_per_subint = pfi.sub.bytes_per_subint * pf.hdr.nbits / pfi.hdr.nbits;
  pf.sub.FITS_typecode = TFLOAT;
  pfi.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
  pfi.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
  pfi.sub.dat_offsets = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
  pfi.sub.dat_scales  = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
  pfi.sub.rawdata = (unsigned char *)malloc(pfi.sub.bytes_per_subint);
  pfi.sub.data = (unsigned char *)malloc(pfi.sub.bytes_per_subint);

#ifdef HAVE_PSRCHIVE  
  // If CAL, construct
  // the archive from which a calibrator will be constructed                                                                                         
  Reference::To<Pulsar::Archive> model_arch;

  // default calibrator type                                                                                                                         
  Reference::To<const Pulsar::Calibrator::Type> pcal_type;
  pcal_type = new Pulsar::CalibratorTypes::SingleAxis;

  if (have_cal) {
      cerr << "Loading calibrator from " << cal_file << endl;

      model_arch = Pulsar::Archive::load(model_file);

      if (model_arch->get<PolnCalibratorExtension>())
	  model_calibrator = new Pulsar::PolnCalibrator(model_arch);
      
      //pcal_type = model_calibrator->get_type();
  }
  
#endif

  // Copy header for the output fil
  memcpy(&pf, &pfi, sizeof(pfi));
  psrfits_create(&pf);


  for (int iloop=0; iloop<pf.hdr.nchan; iloop+=nchanperloop) {
      

      int rv = 0;
      while (run) {

	  // Read Nchans from data
	  rv = psrfits_read_subint(&pfi);

	  if (rv) {
	      run=0;
	      break;
	  }

	  // realloc, copy data

      }

      // Create threads
      for (i=0; i<numfiles; i++) {
	  status = pthread_create(&threads[i], NULL, &merge_psrfits_thread, &fargs[i]);
      }


      // Waiting for threads to finish
      for (i=0; i<numfiles; i++) {
	  pthread_join(threads[i], NULL);
	  if (cmd->verboseP) {
	      for(j=0;j<fargs[i].pf.hdr.nchan; j++) printf("File #%d  Freq[%03d]=%f\n", i, j, fargs[i].pf.sub.dat_freqs[j]);
	  }
      }
      for (i=0; i<numfiles; i++)  statsum += fargs[i].status;
      if (statsum) break;

  }


#if 0
	
    // This is what you would update for each time sample (likely just
    // adjusting the pointer to point to your data)
    pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);
	
    // Here is the real data-writing loop
    do {
      // Update the pf.sub entries here for each subint
      // as well as the pf.sub.data pointer
      nspec = fread(pf.sub.rawdata, pf.hdr.nbits * pf.hdr.nchan * pf.hdr.npol/8, pf.hdr.nsblk, pfi);

		
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

    printf("Done.  Wrote %d subints (%f sec) in %d files.  status = %d\n", 
           pf.tot_rows, pf.T, pf.filenum, pf.status);

#endif

    exit(0);
}
