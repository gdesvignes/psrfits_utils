#include "psrfits.h"
typedef struct {
    //int status;
    //int chan_id;
    //char filename[128];
    int clean;
    unsigned int numdata;
    float *data;


} thread_args;

void *calclean_psrfits_thread(void *args);

pthread_mutex_t lock_calclean;
