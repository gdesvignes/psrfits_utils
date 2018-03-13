double delay_from_dm(double dm, double freq_emitted);

void dedisperse(float* inbuffer,
                float* outbuffer,
                int* delays,
                int maxdelay,
                int nchans,
                int nsamps,
                int index );
