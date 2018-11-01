CC=cc
LAPACK=-llapacke -lblas -lm
TEST=test.o tools.o initial.o
MODULE=module.o tools.o initial.o core.o
test : $(TEST)
	cc -o a.out $(TEST) $(LAPACK) 
mtest : $(MODULE)
	cc -g -o a.out $(MODULE) $(LAPACK) 
initial.o : initial.c
	cc -c initial.c -o initial.o
tools.o : tools.c
	cc -c tools.c -o tools.o
test.o : test.c
	cc -c test.c -o test.o
module.o : module_test.c
	cc -c module_test.c -o module.o
core.o : core.c
	cc -c core.c -o core.o $(LAPACK)
clean :
	rm a.out *.o
