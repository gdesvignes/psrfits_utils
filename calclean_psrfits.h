//#ifndef CALCLEAN_PSRFITS_H
//#define CALCLEAN_PSRFITS_H
#include <vector>
#include <Jones.h>
#include "psrfits.h"
#include "ekstrom.h"

typedef struct {
  int ithread;
  int clean_baseline;
  int redn_clean;
  int median_clean;
  unsigned long totnpts;
  int nchan_per_thread;
  int npol;
  int nchan;
  int nbits;
  int nskip;
  unsigned char *pfiraw;
  unsigned char *pfraw;
  RunningMean<float> *rm;
  std::vector< Jones<float> > response;

} thread_args;


void *calclean_psrfits_32b_thread(void *args);

void *calclean_psrfits_8b_thread(void *args);
//pthread_mutex_t lock_calclean;

//#endif
