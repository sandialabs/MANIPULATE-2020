SHELL=/bin/sh

SRCS= filein.f spectra_format.f grspin.f fold.f fitmd.f\
      terp1.f manipulate.f prune.f grprin.f stat.f covlate.f\
	covread.f	vif.f	force_upper.f	correct.f	spectra_out.f\
	filein_response.f	eigen_out.f	balanc.f	balbak.f \
	cdiv.f	fdump.f	hqr.f	hqr2.f	i1mach.f	j4save.f	\
	orthes.f	ortran.f	scopy.f	scopym.f	sgeev.f	\
	xercnt.f	xerhlt.f	xermsg.f	xerprn.f	xersve.f 	df.f\
	xgetua.f	cov_norm.f	cov_combine.f	pka_grprin.f	integral_pka_grprin.f

OBJS= filein.o spectra_format.o grspin.o fold.o fitmd.o\
      terp1.o manipulate.o prune.o grprin.o stat.o covlate.o\
	covread.o	vif.o	force_upper.o	correct.o	spectra_out.o\
	filein_response.o	eigen_out.o	balanc.o	balbak.o \
	cdiv.o	fdump.o	hqr.o	hqr2.o	i1mach.o	j4save.o	\
	orthes.o	ortran.o	scopy.o	scopym.o	sgeev.o	\
	xercnt.o	xerhlt.o	xermsg.o	xerprn.o	xersve.o  	df.o\
	xgetua.o	cov_norm.o	cov_combine.o	pka_grprin.o	integral_pka_grprin.o

FFT=	f95
CFT=	f95
LDR=	f95
LDFLAGS=-g
LIBS=	
COMMAND=../bin/xmanipulate-2020
PROFLIB=
#
# To compile with FlowTrace turned on Comment out the next line (Put
# a '#' in front of it)
# WARNING:  DO NOT PROFILE AND FLOWTRACE AT THE SAME TIME
#
FFLAGS=-fno-range-check -fdec
#
# To compile with FlowTrace turned on unComment the next line (Take
# out the '#' in front of it) and unComment the line following it
# on Cray-2 systems
# WARNING:  DO NOT PROFILE AND FLOWTRACE AT THE SAME TIME
#
#FFLAGS=-g -i4 -lV77 -Nl30 -xl -e f
#LIBS=	 -lflow
#
# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.
#
$(COMMAND):	$(OBJS)
	$(LDR) $(LDFLAGS) -o $(COMMAND) $(OBJS) $(LIBS)
#
# Make the profiled version of the command and call it a.out.prof
# WARNING:  DO NOT PROFILE AND FLOWTRACE AT THE SAME TIME
#
$(COMMAND).prof:	$(OBJS)
	$(LDR) $(LDFLAGS) -o $(COMMAND) $(OBJS) $(PROFLIB) $(LIBS)
clean:
	rm -f *.o *.l
clobber:	clean
	rm -f $(COMMAND)
void:	clobber
	rm -f $(SRCS) makefile
