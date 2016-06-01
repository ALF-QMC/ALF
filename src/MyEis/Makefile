LIB=libeis.a
OBJS=balanc.o  cdiv.o    comqr2.o  eltran.o  htribk.o  rs.o\
tred2.o balbak.o  cg.o      corth.o   epslon.o  htridi.o  tql2.o\
cbabk2.o  ch.o      csroot.o  hqr.o     pythag.o  tqlrat.o\
cbal.o    comqr.o   elmhes.o  hqr2.o    rg.o      tred1.o
$(LIB): $(OBJS)
	ar r $(LIB) $(OBJS)
.f.o:
	$(FC) $(FLAGS) $< 
clean:
	rm $(LIB) $(OBJS)
