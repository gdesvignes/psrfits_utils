#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fbk_header.h"
#include <fcntl.h>
#include <sys/stat.h>

void angle_split(double angle, int *dd, int *mm, double *ss) {
  int negative;
  if (angle<0.0) {
    angle*=-1.0;
    negative=1;
  } else {
    negative=0;
  }
  *dd=(int) (angle/10000.0);
  angle-=(double) (*dd)*10000.0;
  *mm=(int) (angle/100.0);
  *ss=angle-100.0*(*mm);
  if (negative) *dd = *dd * -1;
}

long long sizeof_file(char name[]) /* includefile */
{
  struct stat stbuf;

  if(stat(name,&stbuf) == -1)
    {
      fprintf(stderr, "f_siz: can't access %s\n",name);
      exit(0);
    }
  return(stbuf.st_size);
}

long long nsamples(char *filename,int headersize, header *h) /*includefile*/
{
  long long datasize,numsamps;
  datasize=sizeof_file(filename)-headersize;
  numsamps=(long long) (long double) (datasize)/ (((long double) h->nbits) / 8.0)
    /(long double) h->nifs/(long double) h->nchans;
  return(numsamps);
}


int strings_equal (char const *string1, char const *string2) {
  if (!strcmp(string1,string2)) {
    return 1;
  } else {
    return 0;
  }
}

void get_string(FILE *inputfile, int *nbytes, char string[]) {
  int nchar;
  strcpy(string,"ERROR");
  fread(&nchar, sizeof(int), 1, inputfile);
  if (feof(inputfile)) exit(0);
  if (nchar>80 || nchar<1) return;
  *nbytes=sizeof(int);
  fread(string, nchar, 1, inputfile);
  string[nchar]='\0';
  *nbytes+=nchar;
}



int read_header(FILE *inputfile, header *h) {
  char string[80], message[80];
  int itmp,nbytes,totalbytes,expecting_rawdatafile=0,expecting_source_name=0; 
  int expecting_frequency_table=0,channel_index;


  /* try to read in the first line of the header */
  get_string(inputfile,&nbytes,string);
  if (!strings_equal(string,"HEADER_START")) {
	/* the data file is not in standard format, rewind and return */
	rewind(inputfile);
	return 0;
  }
  /* store total number of bytes read so far */
  totalbytes=nbytes;

  /* loop over and read remaining header lines until HEADER_END reached */
  while (1) {
    get_string(inputfile,&nbytes,string);
    if (strings_equal(string,"HEADER_END")) break;
    totalbytes+=nbytes;
    if (strings_equal(string,"rawdatafile")) {
      expecting_rawdatafile=1;
    } else if (strings_equal(string,"source_name")) {
      expecting_source_name=1;
    } else if (strings_equal(string,"FREQUENCY_START")) {
      expecting_frequency_table=1;
      channel_index=0;
    } else if (strings_equal(string,"FREQUENCY_END")) {
      expecting_frequency_table=0;
    } else if (strings_equal(string,"az_start")) {
      fread(&h->az_start,sizeof(h->az_start),1,inputfile);
      totalbytes+=sizeof(h->az_start);
    } else if (strings_equal(string,"za_start")) {
      fread(&h->za_start,sizeof(h->za_start),1,inputfile);
      totalbytes+=sizeof(h->za_start);
    } else if (strings_equal(string,"src_raj")) {
      fread(&h->src_raj,sizeof(h->src_raj),1,inputfile);
      totalbytes+=sizeof(h->src_raj);
    } else if (strings_equal(string,"src_dej")) {
      fread(&h->src_dej,sizeof(h->src_dej),1,inputfile);
      totalbytes+=sizeof(h->src_dej);
    } else if (strings_equal(string,"tstart")) {
      fread(&h->tstart,sizeof(h->tstart),1,inputfile);
      totalbytes+=sizeof(h->tstart);
    } else if (strings_equal(string,"tsamp")) {
      fread(&h->tsamp,sizeof(h->tsamp),1,inputfile);
      totalbytes+=sizeof(h->tsamp);
    } else if (strings_equal(string,"period")) {
      fread(&h->period,sizeof(h->period),1,inputfile);
      totalbytes+=sizeof(h->period);
    } else if (strings_equal(string,"fch1")) {
      fread(&h->fch1,sizeof(h->fch1),1,inputfile);
      totalbytes+=sizeof(h->fch1);
    } else if (strings_equal(string,"foff")) {
      fread(&h->foff,sizeof(h->foff),1,inputfile);
      totalbytes+=sizeof(h->foff);
    } else if (strings_equal(string,"nchans")) {
      fread(&h->nchans,sizeof(h->nchans),1,inputfile);
      totalbytes+=sizeof(h->nchans);
    } else if (strings_equal(string,"telescope_id")) {
      fread(&h->telescope_id,sizeof(h->telescope_id),1,inputfile);
      totalbytes+=sizeof(h->telescope_id);
    } else if (strings_equal(string,"machine_id")) {
      fread(&h->machine_id,sizeof(h->machine_id),1,inputfile);
      totalbytes+=sizeof(h->machine_id);
    } else if (strings_equal(string,"data_type")) {
      fread(&h->data_type,sizeof(h->data_type),1,inputfile);
      totalbytes+=sizeof(h->data_type);
    } else if (strings_equal(string,"ibeam")) {
      fread(&h->ibeam,sizeof(h->ibeam),1,inputfile);
      totalbytes+=sizeof(h->ibeam);
    } else if (strings_equal(string,"nbeams")) {
      fread(&h->nbeams,sizeof(h->nbeams),1,inputfile);
      totalbytes+=sizeof(h->nbeams);
    } else if (strings_equal(string,"nbits")) {
      fread(&h->nbits,sizeof(h->nbits),1,inputfile);
      totalbytes+=sizeof(h->nbits);
    } else if (strings_equal(string,"barycentric")) {
      fread(&h->barycentric,sizeof(h->barycentric),1,inputfile);
      totalbytes+=sizeof(h->barycentric);
    } else if (strings_equal(string,"pulsarcentric")) {
      fread(&h->pulsarcentric,sizeof(h->pulsarcentric),1,inputfile);
      totalbytes+=sizeof(h->pulsarcentric);
    } else if (strings_equal(string,"nbins")) {
      fread(&h->nbins,sizeof(h->nbins),1,inputfile);
      totalbytes+=sizeof(h->nbins);
    } else if (strings_equal(string,"nsamples")) {
      /* read this one only for backwards compatibility */
      fread(&itmp,sizeof(itmp),1,inputfile);
      totalbytes+=sizeof(itmp);
    } else if (strings_equal(string,"nifs")) {
      fread(&h->nifs,sizeof(h->nifs),1,inputfile);
      totalbytes+=sizeof(h->nifs);
    } else if (strings_equal(string,"npuls")) {
      fread(&h->npuls,sizeof(h->npuls),1,inputfile);
      totalbytes+=sizeof(h->npuls);
    } else if (strings_equal(string,"refdm")) {
      fread(&h->refdm,sizeof(h->refdm),1,inputfile);
      totalbytes+=sizeof(h->refdm);
    } else if (expecting_rawdatafile) {
      strcpy(h->rawdatafile,string);
      expecting_rawdatafile=0;
    } else if (expecting_source_name) {
      strcpy(h->source_name,string);
      expecting_source_name=0;
    } else if (strings_equal(string,"smin")) {
      fread(&h->smin,sizeof(h->smin),1,inputfile);
      totalbytes+=sizeof(h->smin);
    } else if (strings_equal(string,"smax")) {
      fread(&h->smax,sizeof(h->smax),1,inputfile);
      totalbytes+=sizeof(h->smax);
    } else {
      sprintf(message,"read_header - unknown parameter: %s\n",string);
      fprintf(stderr,"ERROR: %s\n",message);
      exit(1);
    } 
  } 

  /* add on last header string */
  totalbytes+=nbytes;

  /* return total number of bytes read */
  return totalbytes;
}
