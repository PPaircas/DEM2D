# compile library file 
CFILE1=basefun
LIBFILE=libmyfun

# Compile path
GCC = /usr/bin/gcc


all:
	-$(GCC) -c $(CFILE1).c
	ar rcs $(LIBFILE).a $(CFILE1).o
	rm -f *.o

clean::
	rm -f $(LIBFILE).a
