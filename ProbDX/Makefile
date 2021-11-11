CC = icpc
CFLAGS = -O3 -std=c++0x -I/usr/include/x86_64-linux-gnu/c++/4.8/ -openmp -I../GeomBD/
#CFLAGS += -DDEBUG
#CFLAGS += -DTOUCHING
LDFLAGS = -lcilkrts -openmp
		
O_gridder = Main.o 
X_gridder = ProbDX
		

all : gridder 

gridder: $(O_gridder) $(C_gridder)
	$(CC) $(CFLAGS) $(O_gridder) $(C_gridder) -o $(X_gridder) $(LDFLAGS)


%.o : %.cc
	@echo " - Compiling $<..." 
	$(CC) $(CFLAGS) -c $<

%.o : %.cu
	@echo " - Compiling $<..." 
	$(NVCC) -c $<

clean:
	@rm *.o $(X_gridder) $(X_griddes) $(X_bpm2dx) $(X_griddex)
