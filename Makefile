headers = ../Include/
                                                                                
#CC = /usr/apps/gcc/4.3.2/bin/g++
CC = /usr/bin/g++ -g -Wall 
#-O2 -funroll-loops -ffast-math -fomit-frame-pointer
               

AR = ar ru
                                                                 
LIB = nldr_lib.a
LIBDIR = ../../Libraries


OBJS = nldr_Isomap.o \
	nldr_kpca.o \
	nldr_LapEigMap.o \
	nldr_LLEmbed.o \
	nldr_LTSA.o \
	nldr_NLDimRed.o \
	nldr_pca.o \
	nldr_tSNE.o \
	
	

library: libxx
	$(AR) $(LIB) $(OBJS)
	ranlib $(LIB)
	cp $(LIB) $(LIBDIR)

libxx: $(OBJS)



SUFFIXES = .cxx .C
.SUFFIXES: $(SUFFIXES)
                                                                                
.cxx.o:
	$(CC) -c -I. -I$(headers) $<
.cc.o:	
	$(CC) -c $<

.C.o:	
	$(CC) -c $<
clean:
	$(RM) -R *.o *.a core 
