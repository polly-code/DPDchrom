objects = dpdnanov1.25.F90

   FC = mpif90		
   FFLAGS = 
   LDFLAGS = -O3 -ffree-line-length-none
   TARGET = dpd

default: $(objects) 
	$(FC) $(LDFLAGS) -o $(TARGET) $(objects)
   $(objects) :

clean: 
	rm -f $(TARGET) *.o
