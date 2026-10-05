include ../SLmake.inc

TESTOBJS = psscaex.o pdscaexinfo.o
TESTOBJD = pdscaex.o pdscaexinfo.o
TESTOBJC = pcscaex.o pdscaexinfo.o
TESTOBJZ = pzscaex.o pdscaexinfo.o
TESTOBJG = pdgaussex.o

all: pdscaex pcscaex pzscaex psscaex pdgaussex

single: psscaex

double: pdscaex pdgaussex

complex: pcscaex

complex16: pzscaex

psscaex: $(TESTOBJS)
	$(FCLOADER) $(FCLOADFLAGS) -o xsscaex $(TESTOBJS) ../$(SCALAPACKLIB) $(LIBS)

pzscaex: $(TESTOBJZ)
	$(FCLOADER) $(FCLOADFLAGS) -o xzscaex $(TESTOBJZ) ../$(SCALAPACKLIB) $(LIBS)

pcscaex: $(TESTOBJC)
	$(FCLOADER) $(FCLOADFLAGS) -o xcscaex $(TESTOBJC) ../$(SCALAPACKLIB) $(LIBS)

pdscaex: $(TESTOBJD)
	$(FCLOADER) $(FCLOADFLAGS) -o xdscaex $(TESTOBJD) ../$(SCALAPACKLIB) $(LIBS)

pdgaussex: $(TESTOBJG)
	$(FCLOADER) $(FCLOADFLAGS) -o xdgaussex $(TESTOBJG) ../$(SCALAPACKLIB) $(LIBS)

clean :
	rm -f $(TESTOBJS) $(TESTOBJD) $(TESTOBJZ) $(TESTOBJC) $(TESTOBJG) xsscaex xzscaex xcscaex xdscaex \
	      xdgaussex

.f.o : ; $(FC) -c $(FCFLAGS) $*.f

.c.o : ; $(CC) -c $(CDEFS) $(CCFLAGS) $*.c

