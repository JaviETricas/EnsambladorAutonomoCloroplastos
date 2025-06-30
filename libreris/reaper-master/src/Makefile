



# NOTE.
# SWNUM needs to be compiled into various programs. Caveat.



CC = gcc
COPTIM = -Wall -O3 -std=gnu99 -DHAVE_ZLIB=1 -DTRACING_ON=0
	#
	# -march=native -mtune=native
	# icc options for optimisation.  COPTIM = -xHOST -O3 -no-prec-div -static -DHAVE_ZLIB=1

CDEBUG = -g -Wall
CARRRR = -Wall -O3 -fPIC -std=gnu99
CBIG   = -Wall -O2 -g -DMAXFIELDSIZE=255
CSTATIC= -Wall -O3 -std=gnu99 -DHAVE_ZLIB=1 -static
CSMALL = -Wall -Os -std=gnu99 -DHAVE_ZLIB=1 -fpack-struct


CMAIN=$(CDEBUG)
L = -lz -lm
RI = -I${R_INCLUDE_DIR} -I${R_INCLUDE_DIR}/${R_ARCH}
RI = -I${R_INCLUDE_DIR}


ifndef TSW
	TSW = double
endif

ifndef T
	T = optim
endif


ifdef R_INCLUDE_DIR
	T = reaper.so
	.DEFAULT_GOAL := reaper.so
endif


SKT = unsigned long
ifeq ("$(SWANINDEXBITS)", "128")
   SKT = unsigned __int128
endif


ifeq ("$(T)", "optim")
	CMAIN = $(COPTIM)
else ifeq ("$(T)", "big")
	CMAIN = $(CBIG)
else ifeq ("$(T)", "lean")
	CMAIN = $(CSMALL)
else ifeq ("$(T)", "debug")
	CMAIN = $(CDEBUG)
else ifeq ("$(T)", "static")
	CMAIN = $(CSTATIC)
else ifeq ("$(T)", "reaper.so")
	CMAIN = $(CARRRR)
else
	CMAIN = $(COPTIM)
endif


ifeq ("$(TSW)", "unsigned")
	sw_format = unsigned
else ifeq ("$(TSW)", "char")
	sw_format = "unsigned char"
else ifeq ("$(TSW)", "double")
	sw_format = double
else ifeq ("$(TSW)", "int")
	sw_format = int
else ifeq ("$(TSW)", "")
	sw_format = double
else
   $(error Not supported type for TSW: $(TSW))
endif

CSW = -DSWNUM="$(sw_format)" -DSWAN_K_TYPE="$(SKT)"


all : reaper tally minion swan transpose
.PHONY : all


reaper: reaper.o sw.o slib.o table.o trint.o version.o dna.o
	$(CC) -fstrict-aliasing $(CMAIN) -o reaper trint.o version.o dna.o reaper.o sw.o table.o slib.o $(L)

swan: swan.o sw.o slib.o dna.o version.o
	$(CC) -fstrict-aliasing $(CMAIN) -o swan slib.o swan.o sw.o dna.o version.o $(L)

transpose: transpose.o sw.o slib.o dna.o version.o table.o
	$(CC) -fstrict-aliasing $(CMAIN) -o transpose slib.o transpose.o sw.o dna.o version.o table.o $(L)

sword: sword.o sw.o slib.o dna.o version.o table.o
	$(CC) -fstrict-aliasing $(CMAIN) -o sword slib.o sword.o sw.o dna.o version.o table.o $(L)

tally: tally.o slib.o trint.o version.o dna.o sw.o
	$(CC) -fstrict-aliasing $(CMAIN) -o tally trint.o dna.o tally.o version.o slib.o $(L)

minion: minion.o slib.o version.o trint.o sw.o
	$(CC) -fstrict-aliasing $(CMAIN) -o minion version.o trint.o slib.o minion.o sw.o $(L)

reaper.so: .clean bindings.o tally-R.o reaper-R.o sw.o slib.o table.o trint.o version.o
	R CMD SHLIB -o reaper.so bindings.o reaper-R.o tally-R.o dna.o sw.o slib.o table.o trint.o version.o
	# $(CC) -std=gnu99 $(CMAIN) -shared -o reaper.so bindings.o sw.o slib.o table.o trint.o


tmptst:
	$(CC) -DRREAPERR -fstrict-aliasing -Wall -O2 -fPIC -c -o bindings.o reaper.c -I/nfs/research2/enright/local/lib64/R/include -I/nfs/research2/enright/local/lib64/R/include/

mytest:
	echo .target=$(.TARGET)
	echo CMAIN=$(CMAIN)
	echo T=$(T)
	echo TSW=$(TSW) CSW=$(CSW)


bindings.o: bindings.c reaper.h tally.h
	$(CC) -fstrict-aliasing $(CMAIN) -c -o bindings.o bindings.c $(RI)

chair: chair.o table.o slib.o
	$(CC) -g -fstrict-aliasing $(CMAIN) -o chair chair.o table.o slib.o $(L)

trint.o: trint.c trint.h
	$(CC) -fstrict-aliasing $(CMAIN) -c -o trint.o trint.c

sw.o: sw.c sw.h
	$(CC) -fstrict-aliasing $(CMAIN) $(CSW) -c -o sw.o sw.c

version.o: version.c version.h
	$(CC) -fstrict-aliasing $(CMAIN) -c -o version.o version.c

slib.o: slib.c slib.h
	$(CC) -fstrict-aliasing $(CMAIN) -c -o slib.o slib.c

minion.o: minion.c
	$(CC) -fstrict-aliasing $(CMAIN) $(CSW) -c -o minion.o minion.c

swan.o: swan.c
	$(CC) -fstrict-aliasing $(CMAIN) $(CSW) -c -o swan.o swan.c

transpose.o: transpose.c
	$(CC) -fstrict-aliasing $(CMAIN) $(CSW) -c -o transpose.o transpose.c

sword.o: sword.c
	$(CC) -fstrict-aliasing $(CMAIN) $(CSW) -c -o sword.o sword.c

dna.o: dna.c dna.h
	$(CC) -g -fstrict-aliasing $(CMAIN) -c -o dna.o dna.c

table.o: table.c table.h slib.h
	$(CC) -g -fstrict-aliasing $(CMAIN) -c -o table.o table.c

#reaper.o: reaper.c slib.h sw.c slib.c trint.c table.c
#	$(CC) -fstrict-aliasing $(CMAIN) -c -o reaper.o reaper.c

reaper.o: reaper.c sw.h slib.h trint.h table.h
	$(CC) -fstrict-aliasing $(CMAIN) $(CSW) -c -o reaper.o reaper.c

reaper-R.o: reaper.c sw.h slib.h trint.h table.h
	$(CC) -DBUILD_R_BINDINGS -fstrict-aliasing $(CMAIN) -c -o reaper-R.o reaper.c


tally.o: tally.c
	$(CC) -fstrict-aliasing $(CMAIN) -c -o tally.o tally.c

matecare4.o: matecare4.c
	$(CC) -fstrict-aliasing $(CMAIN) $(CSW) -c -o matecare4.o matecare4.c

tally-R.o: tally.c
	$(CC) -DBUILD_R_BINDINGS -fstrict-aliasing $(CMAIN) -c -o tally-R.o tally.c


wc:
	wc trint.c reaper.c sw.c table.c slib.c tally.c swan.c

clean .clean:
	rm -f {bindings,trint,reaper,reaper-R,slib,sw,table,tally,tally-R,version,minion,dna,swan,sword,transpose.o}.o



matecare4: matecare4.o slib.o trint.o version.o dna.o sw.o
	$(CC) -fstrict-aliasing $(CMAIN) -o matecare4 sw.o trint.o dna.o matecare4.o version.o slib.o $(L)
