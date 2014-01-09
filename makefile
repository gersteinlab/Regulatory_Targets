# Project: DRM2Target
#DRM2Target.cpp tabixSearch.cpp gff.cpp vmstat.cpp
CC   = g++ -g -pg  -fpermissive -std=c++0x  -DVERSION=1.07
INCS = -Ialglib  -Itabix-0.2.6 
LIBS =  -Ltabix-0.2.6
OBJ  = DRM2Target.o tabixSearch.o gff.o vmstat.o
LINKOBJ  = DRM2Target.o tabixSearch.o gff.o vmstat.o
BUILD_IDX_OBJ =
BIN  = DRM2Target
CXXFLAGS = $(CXXINCS)   -fexpensive-optimizations -O3
CFLAGS = $(INCS)   -fexpensive-optimizations -O3
RM = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before DRM2Target all-after


clean: 
	rm -f  *.o  DRM2Target

target: DRM2Target


DRM2Target: $(OBJ)
	$(CC) $(INCS) $(LIBS) $(LINKOBJ) -o $@ -ltabix -lm -lz -lgsl -lgslcblas -lboost_regex alglib.a

DRM2Target.o:DRM2Target.cpp tabixSearch.cpp vmstat.cpp gff.cpp drmPair.h stat.h adjustp.h
	$(CC) $(INCS) $(LIBS) -c DRM2Target.cpp -ltabix -lm -lz -lgsl -lgslcblas -lboost_regex 

vmstat.o: vmstat.cpp
	$(CC) -c vmstat.cpp

gff.o: gff.cpp segment.h
	$(CC) $(LIBS) -c gff.cpp  -lboost_regex

tabixSearch.o: tabixSearch.cpp drmPair.h
	$(CC) $(INCS) $(LIBS) -c tabixSearch.cpp  -ltabix



