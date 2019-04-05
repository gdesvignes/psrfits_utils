#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <getopt.h>
#include <stdbool.h>
#include <inttypes.h>
#include "psrfits.h"
#include "calclean_psrfits.h"

#include <string>
#include <iostream>
#include "ekstrom.h"

// Record the execution time of some code, in milliseconds.
#if 0
#include <time.h>
#define DECLARE_TIMING(s)  clock_t timeStart_##s; double timeDiff_##s; double timeTally_##s = 0; int countTally_##s = 0
#define START_TIMING(s)    timeStart_##s = clock()
#define STOP_TIMING(s)     timeDiff_##s = (double)(clock() - timeStart_##s); timeTally_##s += timeDiff_##s; countTally_##s++
#define GET_TIMING(s)      (double)(timeDiff_##s / (CLOCKS_PER_SEC*1000.0))
#define GET_AVERAGE_TIMING(s)   (double)(countTally_##s ? timeTally_##s/ ((double)countTally_##s * CLOCKS_PER_SEC*1000.0) : 0)
#define CLEAR_AVERAGE_TIMING(s) timeTally_##s = 0; countTally_##s = 0
#endif

#ifdef HAVE_PSRCHIVE
#include "Pulsar/Archive.h"
#include "Pulsar/PolnCalibrator.h"
#include "Pulsar/CalibratorStokes.h"
#include "Pulsar/CalibratorTypes.h"
#include "Pulsar/PolnCalibratorExtension.h"
#endif

//using namespace Pulsar;

void usage() {
  printf(
		 "Usage: dada2psrfits [options] input_filename_base\n"
		 "Options:\n"
		 "  -h, --help               Print this\n"
#ifdef HAVE_PSRCHIVE
		 "  -c, --cal CAL_archive    Calibrate using a PSRchive noise diode observation\n"
#endif
		 "  -m, --median_clean       Clean the baseline using running median of width 'w'\n"
		 "  -n, --ncploop            Process only ncploop channels per loop\n"
		 "  -s, --skip               Skip nsamples for the median\n"
		 "  -t, --nthreads           Use nthreads for processing\n"
		 "  -w, --window             Set the window to compute the running mean\n"
		 );
}

int main(int argc, char *argv[]) {

  /* Cmd line */
  static struct option long_opts[] = {
	{"help",      0, NULL, 'h'},
	{"cal",       1, NULL, 'c'},
	{"median_clean", 0, NULL, 'm'},
	{"ncploop",   1, NULL, 'n'},
	{"nthreads",  1, NULL, 't'},
	{"window",    1, NULL, 'w'},
	{0,0,0,0}
  };

  int status,opt, opti;
  int rv=0,nchan_per_loop=1, nthreads=1;
  int boxsize = 4096, nskip=1, fnum_start=1;
  struct psrfits pfi, pf;
  char ra_str[16], dec_str[16], source[24];
  char cal_file[128];
  bool median_clean=false, have_cal_file=false;
  while ((opt=getopt_long(argc,argv,"c:mn:s:t:w:h", long_opts,&opti))!=-1) {
      switch (opt) {
#ifdef HAVE_PSRCHIVE
      case 'c':
		strncpy(cal_file, optarg, 128);
		have_cal_file = true;
		break;
#endif
	  case 'm':
		median_clean = true;
		break;
      case 'n':
		nchan_per_loop = atoi(optarg);
		break;
	  case 's':
		nskip = atoi(optarg);
		break;
	  case 't':
		nthreads = atoi(optarg);
		break;
	  case 'w':
		boxsize = atoi(optarg);
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

  // Init open file
  pfi.filenum = 1;
  strcpy(pfi.basefilename, argv[optind]);
  pfi.numfiles = 0;
  pfi.status = 0;

  psrfits_set_files(&pfi, argc - optind, argv + optind);

  if (pfi.numfiles==0) pf.filenum = fnum_start;
  pf.tot_rows = pf.N = pf.T = pf.status = 0;
  
  // Open input file
  status = psrfits_open(&pfi);
  if (status) fits_report_error(stderr, status);

  // Copy header for the output file
  memcpy(&pf, &pfi, sizeof(pfi));

  pf.filenum=0;
  // Change output to 32 bits
  if (have_cal_file) {
	pf.hdr.nbits = 32;
	pf.sub.FITS_typecode = TFLOAT;
	std::cout << "Want to calibrate. Use 32 bits "<< std::endl;
  }
  else {
	if (pfi.hdr.nbits == 32) {
	  pf.hdr.nbits = 32;
	  pf.sub.FITS_typecode = TFLOAT;
	}
	else {
	  pf.hdr.nbits = pfi.hdr.nbits;
	  pf.sub.FITS_typecode = TBYTE;
	}
	//std::cout << "Use 4 bits anyway"<< std::endl;
  }
  pf.rows_per_file = INT_MAX;

  // Change output filename
  if (have_cal_file) 
	sprintf(pf.basefilename, "%sCLEANCAL_%s_%05d_%4d", pf.hdr.backend, pf.hdr.source, (int) pf.hdr.MJD_epoch, (int) pf.hdr.fctr);
  else
	sprintf(pf.basefilename, "%sCLEAN_%s_%05d_%4d", pf.hdr.backend, pf.hdr.source, (int) pf.hdr.MJD_epoch, (int) pf.hdr.fctr);

  // Alloc data buffers for the input & output PSRFITS files
  pf.sub.bytes_per_subint = pfi.sub.bytes_per_subint * pf.hdr.nbits / pfi.hdr.nbits;
  pfi.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
  pfi.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
  pfi.sub.dat_offsets = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
  pfi.sub.dat_scales  = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);

  pfi.sub.rawdata = (unsigned char *)malloc(pfi.sub.bytes_per_subint);
  pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);
  //float *pfraw = (float *) pf.sub.rawdata;  
  //unsigned char *pfiraw = (unsigned char *) pfi.sub.rawdata;

  // Default is to process all channels at once
  if (nchan_per_loop==1)
       nchan_per_loop = pf.hdr.nchan;

#ifdef HAVE_PSRCHIVE  
  // If CAL, construct
  // the archive from which a calibrator will be constructed                                                                                         
  Reference::To<Pulsar::Archive> arch;

  // the calibrator constructed from the specified archive
  Reference::To<Pulsar::PolnCalibrator> calibrator;
  
  // default calibrator type
  Reference::To<const Pulsar::Calibrator::Type> pcal_type;
  ///pcal_type = new Pulsar::CalibratorTypes::SingleAxis;

  if (have_cal_file) {
	std::cerr << "Loading calibrator from " << cal_file << std::endl;

	arch = Pulsar::Archive::load(cal_file);

	if (arch->get_type() != Signal::PolnCal) {
	  std::cerr << "Archive " << cal_file << " is not a Calibrator" << std::endl;
	  return 0;
	}
	
	std::cerr << "Here" << std::endl;
	calibrator = new Pulsar::PolnCalibrator(arch);
	std::cerr << "Nchan" << std::endl;
	std::cout << calibrator->get_Archive()->get_centre_frequency() << std::endl;
	  
	  //std::cerr << "Nchan: " << model_calibrator.get_nchan() << std::endl;
      //pcal_type = model_calibrator->get_type();

	std::string pcal_file;
	pcal_file = calibrator->get_filenames();
	std::cout << "pac: PolnCalibrator constructed from:\n\t" << pcal_file << std::endl;
  }
  
#endif

  // Create threads
  pthread_t threads[nthreads];  // Thread ids
  thread_args fargs[nthreads];  // Arguments passed to the threads
  int nchan_per_thread = nchan_per_loop / nthreads;
  //unsigned long totnpts = (unsigned long)pf.tot_rows * pf.hdr.nsblk;

  // Init struct for running median
  std::cout << "Boxsize="<< boxsize<<" Time="<< boxsize * pf.hdr.dt<<std::endl;
  RunningMean<float> *rm;
  rm = new RunningMean<float>[pf.hdr.nchan*pf.hdr.npol];
  for (int i=0; i < pf.hdr.nchan*pf.hdr.npol; i++)
	rm[i].set_boxsize(boxsize); // Set boxsize for running median
  
  for (int i=0; i<nthreads; i++) {
	printf("Setting up thread #%d Offset\n", i);
	fargs[i].ithread = i;
	fargs[i].median_clean = median_clean;
	fargs[i].totnpts = pf.hdr.nsblk;
	fargs[i].nchan_per_thread = nchan_per_thread;
	fargs[i].nchan = pf.hdr.nchan;
	fargs[i].npol = pf.hdr.npol;
	fargs[i].nbits = pf.hdr.nbits;
	fargs[i].nskip = nskip;
	fargs[i].pfiraw = pfi.sub.rawdata;
	fargs[i].pfraw = pf.sub.rawdata;
	fargs[i].rm = &rm[i*nchan_per_thread*pf.hdr.npol];
	fargs[i].response.resize(nchan_per_thread);
	//for (int j=0; j<nchan_per_thread; j++) //
	//fargs[i].response[j] = calibrator->get_response(i*nchan_per_thread + j);
  }
  

  // Create output file
  psrfits_create(&pf);

  int run=1;
  while (run) {
      // Read all subints
      rv = psrfits_read_subint(&pfi);

      if (rv) {
		run=0;
		break;
      }
      
      // Copy subint 
      copy_subint_params(&pf, &pfi);
	  pf.sub.dat_freqs = pfi.sub.dat_freqs;
	  pf.sub.dat_weights = pfi.sub.dat_weights;
	  pf.sub.dat_offsets = pfi.sub.dat_offsets;
	  pf.sub.dat_scales = pfi.sub.dat_scales;

	  printf("Create threads\n");
	  for (int i=0; i<nthreads; i++) {
		if (pfi.hdr.nbits == 8)
		  status = pthread_create(&threads[i], NULL, &calclean_psrfits_8b_thread, &fargs[i]);
		if (pfi.hdr.nbits == 32)
		  status = pthread_create(&threads[i], NULL, &calclean_psrfits_32b_thread, &fargs[i]);
		
	  }

	  // Waiting for threads to finish
	  for (int i=0; i<nthreads; i++) {
		pthread_join(threads[i], NULL);
	  }

	  //for (int k=0;k<100;k++) printf("%d %f\n", k, pfraw[k]);
	  
      // Write to disk, with new 32 bits output
      status = psrfits_write_subint(&pf);
      if (status) {
	  // Check to see if subint was successfully created
		fprintf(stderr, "Error writing subint\n");
		fits_report_error(stderr, status);
		exit(1);
      }

  }

  // Finished rewriting the file - Close it!
  fits_close_file(pf.fptr, &(pf.status));
  
  //printf("Execution time: %f ms.\n", GET_TIMING(myTimer) );
  //printf("Average time: %f ms per iteration.\n", GET_AVERAGE_TIMING(myTimer) );
  
  // Close the last file and cleanup
  free(pfi.sub.dat_freqs);
  free(pfi.sub.dat_weights);
  free(pfi.sub.dat_offsets);
  free(pfi.sub.dat_scales);
  free(pfi.sub.rawdata);
  free(pf.sub.rawdata);
  
  printf("Done.  Wrote %d subints (%f sec) in %d files.  status = %d\n", 
           pf.tot_rows, pf.T, pf.filenum, pf.status);

  exit(0);
}
