LIB=libmscbdecomp.a
OBJS=vertex.o Exponentials.o colorvertex.o MvG_decomp.o
MODS=

$(LIB): $(OBJS)
	ar -r $(LIB) $(OBJS)

%.o: %.F90
	$(ALF_FC) -c $(ALF_FLAGS_MODULES) $<

clean:
	rm -f $(OBJS) $(MODS) $(LIB)
