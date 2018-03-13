// Taken from sigPyProc by E. Barr (I think)

void dedisperse(float* inbuffer,
                float* outbuffer,
                int* delays,
                int maxdelay,
                int nchans,
                int nsamps,
                int index )
{
  int ii,jj,cnt;
#pragma omp parallel for default(shared) private(ii,jj,cnt) shared(outbuffer,inbuffer)
  for (ii=0;ii<(nsamps-maxdelay);ii++){
    cnt=0;
    for (jj=0;jj<nchans;jj++){
      outbuffer[index+ii] += inbuffer[(ii*nchans)+(delays[jj]*nchans)+jj];
      cnt++;
    }
    outbuffer[index+ii] /= (float)cnt;
  }
}

