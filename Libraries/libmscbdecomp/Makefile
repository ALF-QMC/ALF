LIB=libmscbdecomp.a
OBJS=Node.o vertex.o SingleColExpBase.o ZeroDiagSingleColExp.o HomogeneousSingleColExp.o TraceLessSingleColExp.o GeneralSingleColExp.o Exponentials.o colorvertex.o graphdata.o MvG_decomp.o
MODS=

$(LIB): $(OBJS)
	ar -r $(LIB) $(OBJS)

%.o: %.F90
	$(ALF_FC) -Wall -c $(ALF_FLAGS_MODULES) $<

clean:
	rm -f $(OBJS) $(MODS) $(LIB)
