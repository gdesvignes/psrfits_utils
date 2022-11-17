//#ifndef CALCLEAN_PSRFITS_H
//#define CALCLEAN_PSRFITS_H
#include <vector>
#include "Pauli.h"
#include "psrfits.h"

typedef struct {
    int ithread;
    unsigned long totnpts;
    int nchan_per_thread;
    int nchan;
    int npol;
    int nbits;
    float *weights;
    unsigned char *pfiraw;
    unsigned char *pfraw;
    std::vector<Matrix<4,4,double>> response;
    double rcvr_sa;
} thread_args;


void *calibrate_psrfits_32b_thread(void *args);

//pthread_mutex_t lock_calclean;

//#endif
