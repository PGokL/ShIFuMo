PROG = ShIFuMo.x
SRC = ShIFuMo.cu

CUDA = nvcc
FLAGS =

compile: $(PROG)

$(PROG): $(SRC)
	$(CUDA) $(SRC) -o $@

clean:
	rm -rf *.o *.x
