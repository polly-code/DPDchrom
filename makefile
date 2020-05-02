objects = DPDchrom.F90

   FC = mpif90		
   FFLAGS = 
   LDFLAGS = -O3 -ffree-line-length-none
   #LDFLAGS = -O0 -CB -traceback -g
   TARGET = dpdchrom

default: $(objects) 
	$(FC) $(LDFLAGS) -o $(TARGET) $(objects)
   $(objects) :

clean: 
	rm -f $(TARGET) *.o
