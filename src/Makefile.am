AM_CXXFLAGS = -O3
AM_CFLAGS = -O3
LIBS = -lfftw3f_mpi -lfftw3f -lgsl -lgslcblas -lm -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw -lcfitsio
common = arrayoperations.C arrayoperations.h cosmology.C cosmology.h geometry.h\
	globaltablevars.h math.C memorytracking.C memorytracking.h\
        parallel_io.C parallel_io.h tables.C tables.h chealpix.c chealpix.h

bin_PROGRAMS = delta2zreion ics limber flatskymap allskymap
delta2zreion_SOURCES = allocate_d2z.C delta2zreion_d2z.C proto_d2z.h\
	allvars_d2z.C bubblehistory_d2z.C io_d2z.C proto_d2z.h.bak\
	allvars_d2z.h commandline_d2z.C	main_d2z.C $(common)
ics_SOURCES = allocate_ics.C allvars_ics.C allvars_ics.h\
	commandline_ics.C grf_ics.C io_ics.C main_ics.C proto_ics.h $(common)
limber_SOURCES = allocate_lmb.C allvars_lmb.C allvars_lmb.h\
	commandline_lmb.C delta2velocity_lmb.C io_lmb.C\
	main_lmb.C powerspectrum_lmb.C proto_lmb.h $(common)
flatskymap_SOURCES = allocate_map.C allvars_map.C allvars_map.h\
	commandline_map.C delta2velocity_map.C io_map.C\
	main_map.C makemaps_map.C proto_map.h $(common)
allskymap_SOURCES = allocate_asm.C allvars_asm.C allvars_asm.h \
	commandline_asm.C delta2velocity_asm.C geometry.C io_asm.C \
	main_asm.C  makemaps_asm.C parallel_io_asm.c \
	proto_asm.h $(common)
