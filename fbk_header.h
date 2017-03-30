
typedef struct {
  int telescope_id;
  int machine_id;
  int data_type;
  int nchans;
  int nbits;
  int nifs;
  int ibeam;
  int nbeams;
  int barycentric;
  int pulsarcentric;
  int nbins;
  long int npuls;
  long int nsamp;
  double tstart;
  double tsamp;
  double fch1;
  double foff;
  double refdm;
  double period;
  double src_raj;
  double src_dej;
  double az_start;
  double za_start;
  double smin;
  double smax;
  char rawdatafile[80];
  char source_name[80];
} header;


int read_header(FILE *inputfile, header *h);
int write_header(FILE *pfo, header *h);
long long nsamples(char *filename,int headersize, header *h);
void angle_split(double angle, int *dd, int *mm, double *ss);
