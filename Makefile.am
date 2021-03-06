ACLOCAL_AMFLAGS = -I config

pkginclude_HEADERS = \
	ascii_header.h \
	calclean_psrfits.h \
	dedisperse.h \
	fold.h \
	merge_psrfits.h \
	polyco.h \
	polyco_struct.h \
	psrfits.h 

noinst_HEADERS = 

dist_data_DATA = \
	guppi_PSRFITS_v3.4_fold_template.txt \
	guppi_PSRFITS_v3.4_search_template.txt

lib_LTLIBRARIES = libpsrfits_utils.la
libpsrfits_utils_la_SOURCES = \
	ascii_header.c \
	dedisperse.c \
    downsample.c \
	fold.c \
	merge_psrfits_thread.c \
	misc_utils.c \
	plot_psrfits_thread.c \
    polyco.c \
    read_psrfits.c \
    write_psrfits.c \
	read_fbk.c

LDADD = libpsrfits_utils.la @CFITSIO_LIBS@
AM_CPPFLAGS = -DPSRFITS_TEMPLATE_DIR='"$(PSRFITS_UTILS)"' @CFITSIO_CFLAGS@
AM_CXXFLAGS = -std=gnu++11

bin_PROGRAMS = fold_psrfits psrfits_singlepulse psrfits_subband combine_mocks merge_psrfits fits_extract dada2psrfits fbk2psrfits calclean_psrfits combine_S45dada_psrfits combine_Pico fbkbands2psrfits
 
if HAVE_PGPLOT
bin_PROGRAMS += plot_psrfits 
endif

fold_psrfits_CFLAGS = @PTHREAD_CFLAGS@ $(AM_CFLAGS)
fold_psrfits_LDFLAGS = @PTHREAD_CFLAGS@ $(AM_LDFLAGS)
fold_psrfits_LDADD = libpsrfits_utils.la @PTHREAD_LIBS@ @CFITSIO_LIBS@

psrfits_subband_SOURCES = psrfits_subband.c \
			  psrfits_subband_cmd.h psrfits_subband_cmd.c

combine_mocks_SOURCES = combine_mocks.c \
				   combine_mocks_cmd.h combine_mocks_cmd.c

merge_psrfits_SOURCES = merge_psrfits.c \
                          merge_psrfits_cmd.h merge_psrfits_cmd.c

merge_psrfits_CFLAGS = @PTHREAD_CFLAGS@ $(AM_CFLAGS)
merge_psrfits_LDFLAGS = @PTHREAD_CFLAGS@ $(AM_LDFLAGS)
merge_psrfits_LDADD = libpsrfits_utils.la @PTHREAD_LIBS@ @CFITSIO_LIBS@


calclean_psrfits_SOURCES = calclean_psrfits.C calclean_psrfits_32b_thread.C calclean_psrfits_8b_thread.C

calclean_psrfits_CPPFLAGS = @PTHREAD_CFLAGS@ $(AM_CFLAGS)
calclean_psrfits_LDFLAGS = @PTHREAD_CFLAGS@ $(AM_LDFLAGS)
calclean_psrfits_LDADD = libpsrfits_utils.la @PTHREAD_LIBS@

if HAVE_PSRCHIVE
calclean_psrfits_CPPFLAGS += @PSRCHIVE_CFLAGS@
calclean_psrfits_LDADD += @PSRCHIVE_LIBS@ 
endif

if HAVE_PGPLOT
plot_psrfits_SOURCES = plot_psrfits.c plot_psrfits_cmd.h plot_psrfits_cmd.c
plot_psrfits_CFLAGS = @PTHREAD_CFLAGS@ @PGPLOT_CFLAGS@ $(AM_CFLAGS)
plot_psrfits_LDFLAGS = @PTHREAD_CFLAGS@ $(AM_LDFLAGS)
plot_psrfits_LDADD = libpsrfits_utils.la @PTHREAD_LIBS@ @PGPLOT_LIBS@ @CFITSIO_LIBS@

if HAVE_PRESTO
plot_psrfits_CFLAGS += @PRESTO_CFLAGS@
plot_psrfits_LDADD += @PRESTO_LIBS@
endif
endif

fits_extract_SOURCES = fits_extract.c
