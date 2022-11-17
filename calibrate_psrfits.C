#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <getopt.h>
#include <stdbool.h>
#include <inttypes.h>
#include "psrfits.h"
#include "calibrate_psrfits.h"

#include <string>
#include <iostream>
#include <vector>

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

#include "Pauli.h"
#include "Pulsar/Archive.h"
#include "Pulsar/PolnCalibrator.h"
#include "Pulsar/PolnCalibratorExtension.h"
#include "Pulsar/SingleAxisCalibrator.h"
#include "Pulsar/Telescope.h"
#include "Pulsar/Receiver.h"
#include "Horizon.h"
#include "MJD.h"

void usage() {
  printf(
		 "Usage: calibrate_psrfits [options] input_filename_base\n"
		 "Options:\n"
		 "  -h, --help               Print this\n"
		 "  -c, --cal CAL_archive    Calibrate using a PSRchive noise diode observation\n"
		 "  -t, --nthreads           Use nthreads for processing\n"
		 );
}

int main(int argc, char *argv[]) {

  /* Cmd line */
  static struct option long_opts[] = {
	{"help",      0, NULL, 'h'},
	{"cal",       1, NULL, 'c'},
	{"nthreads",  1, NULL, 't'},
	{0,0,0,0}
  };

  int status,opt, opti;
  int rv=0,nchan_per_loop=1, nthreads=1, fnum_start=1;
  struct psrfits pfi, pf;
  char ra_str[16], dec_str[16], source[24];
  char cal_file[128];
  bool have_cal_file=false;
  double PA;
  long double mjd;

  Matrix<4,4,double> MPA;
  // Form the parallactic angle matrix 
  MPA[0][0] = 1; MPA[3][3] = 1;
  
  while ((opt=getopt_long(argc,argv,"c:t:h", long_opts,&opti))!=-1) {
      switch (opt) {
      case 'c':
	  strncpy(cal_file, optarg, 128);
	  have_cal_file = true;
	  break;
      case 't':
	  nthreads = atoi(optarg);
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

  // MJD
  mjd = pf.hdr.MJD_epoch;
      
  pf.filenum=0;
  // Change output to 32 bits
  if (have_cal_file) {
	pf.hdr.nbits = 32;
	pf.sub.FITS_typecode = TFLOAT;
	std::cout << "Want to calibrate. Use 32 bits "<< std::endl;
	sprintf(pf.basefilename, "%sCAL_%s_%05d_%4d", pf.hdr.backend, pf.hdr.source, (int) pf.hdr.MJD_epoch, (int) pf.hdr.fctr);
  } else {
	sprintf(pf.basefilename, "%sCLEAN_%s_%05d_%4d", pf.hdr.backend, pf.hdr.source, (int) pf.hdr.MJD_epoch, (int) pf.hdr.fctr);
  }
  
  // Alloc data buffers for the input & output PSRFITS files
  pf.sub.bytes_per_subint = pfi.sub.bytes_per_subint * pf.hdr.nbits / pfi.hdr.nbits;
  pfi.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
  pfi.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
  pfi.sub.dat_offsets = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
  pfi.sub.dat_scales  = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
  pfi.sub.rawdata = (unsigned char *)malloc(pfi.sub.bytes_per_subint);
  
  pf.sub.dat_freqs = (float *)malloc(sizeof(float) * pf.hdr.nchan);
  pf.sub.dat_weights = (float *)malloc(sizeof(float) * pf.hdr.nchan);
  pf.sub.dat_offsets = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
  pf.sub.dat_scales  = (float *)malloc(sizeof(float) * pf.hdr.nchan * pf.hdr.npol);
  pf.sub.rawdata = (unsigned char *)malloc(pf.sub.bytes_per_subint);


  Reference::To<Pulsar::Archive> arch;

  // the calibrator constructed from the specified archive
  Reference::To<Pulsar::PolnCalibrator> calibrator;

  std::vector<Matrix<4,4,double>> response;
  std::vector<int> zero_weight_chans;
  std::vector<int> is_chan_valid;
  Horizon horizon;
  double rcvr_sa;

  if (have_cal_file) {
	std::cerr << "Loading calibrator from " << cal_file << std::endl;

	arch = Pulsar::Archive::load(cal_file);

	// Recompute Parallactic angle from the archive
	const Pulsar::Telescope* telescope = arch->get<Pulsar::Telescope>();
	horizon.set_observatory_latitude (telescope->get_latitude().getRadians());
	horizon.set_observatory_longitude (telescope->get_longitude().getRadians());
	horizon.set_source_coordinates(sky_coord(arch->get_coordinates()));

	// Get Symmetry angle
	const Pulsar::Receiver* receiver = arch->get<Pulsar::Receiver>();
	rcvr_sa = receiver->get_orientation().getRadians();
	
	calibrator = new Pulsar::SingleAxisCalibrator (arch);
	
	for (unsigned ichan=0; ichan<calibrator->get_nchan(); ichan++) {
	    Matrix<4,4,double> M;
	    if (calibrator->get_transformation_valid (ichan)) {
		const MEAL::Complex2* xform = calibrator->get_transformation(ichan);
		Jones<double> J = xform->evaluate();

		M = Mueller (J);
		is_chan_valid.push_back(1);
	    }
	    // Flag channels with invalid calibration
	    else {
		//pf.sub.dat_weights[ichan] = 0.0;
		zero_weight_chans.push_back(ichan);
		is_chan_valid.push_back(0);
	    }
	    response.push_back(M);
	}
  }
  
  // Create threads
  pthread_t threads[nthreads];  // Thread ids
  thread_args fargs[nthreads];  // Arguments passed to the threads
  int nchan_per_thread = pf.hdr.nchan / nthreads;
  //unsigned long totnpts = (unsigned long)pf.tot_rows * pf.hdr.nsblk;

  for (int i=0; i<nthreads; i++) {
	printf("Setting up thread #%d Offset\n", i);
	fargs[i].ithread = i;
	fargs[i].totnpts = pf.hdr.nsblk;
	fargs[i].nchan_per_thread = nchan_per_thread;
	fargs[i].nchan = pf.hdr.nchan;
	fargs[i].npol = pf.hdr.npol;
	fargs[i].nbits = pf.hdr.nbits;
	fargs[i].pfiraw = pfi.sub.rawdata;
	fargs[i].pfraw = pf.sub.rawdata;
	fargs[i].response = response;
	fargs[i].weights = pf.sub.dat_weights;
	fargs[i].rcvr_sa = rcvr_sa;
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
      for (int ichan : zero_weight_chans) pf.sub.dat_weights[ichan] = 0.0;

      // MJD
      mjd += 0.5*pf.sub.tsubint/86400.;
      horizon.set_epoch (MJD((double)mjd));
      PA = horizon.get_parallactic_angle ();

      // Forms the Parallactic angle matrix
      MPA[1][1] = cos(2*PA); MPA[2][2] = cos(2*PA);
      MPA[1][2] = sin(2*PA); MPA[2][1] = -sin(2*PA);

      for (int ithread=0; ithread<nthreads; ithread++) {
	  for (int i=0; i<pf.hdr.nchan; i++) {
	      if (is_chan_valid[i] == 1)
	      fargs[ithread].response[i] = inv(response[i] * MPA);
	  }
      }
      
      //printf("Create threads\n");
      for (int i=0; i<nthreads; i++) {
	  status = pthread_create(&threads[i], NULL, &calibrate_psrfits_32b_thread, &fargs[i]);
	  
      }
      
      // Waiting for threads to finish
      for (int i=0; i<nthreads; i++) {
	  pthread_join(threads[i], NULL);
      }
      
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
