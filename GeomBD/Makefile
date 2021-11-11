CC = icpc
CFLAGS = -O3 -std=c++0x -I/usr/include/x86_64-linux-gnu/c++/4.8/
CFLAGS += -DDEBUG
#CFLAGS += -DTOUCHING
LDFLAGS = -lcilkrts -lmkl_blas95_ilp64 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lmkl_gf_ilp64 -openmp
		
O_bd = Main.o Model.o Model_Input.o Model_Force.o Model_Output.o Body.o Session.o
X_bd = ../bin/GBD

all : BD

BD: $(O_bd) $(C_bd)
	$(CC) $(CFLAGS) $(O_bd) $(C_bd) -o $(X_bd) $(LDFLAGS)

%.o : %.cc
	@echo " - Compiling $<..." 
	$(CC) $(CFLAGS) -c $<

%.o : %.cu
	@echo " - Compiling $<..." 
	$(NVCC) -c $<

clean:
	@rm *.o $(X_bd)
