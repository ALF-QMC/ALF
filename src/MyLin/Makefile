LIB=liblin.a
OBJS= cgedi.o   cgefa.o  dgedi.o dgefa.o zgedi.o zgefa.o zqrdc.o zqrsl.o
$(LIB): $(OBJS)
	ar r $(LIB) $(OBJS)
.f.o:
	$(FC) $(FLAGS) $< 
clean:
	rm $(LIB) $(OBJS)
