CC = g++

VC_TOP     = ../../../../vcollide-2.01

CFLAGS = -w -O

LDFLAGS   = -L$(VC_TOP)/lib -L$(VC_TOP)/RAPID

INCLUDES  = -I$(VC_TOP)/include

all : pearMix_hard.cpp
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o pearMix_hard pearMix_hard.cpp -lVCollide -lRAPID -lm -DPOLY

MONO : pearMix_hard.cpp
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o pearMix_hard pearMix_hard.cpp -lVCollide -lRAPID -lm

clean: 
	rm -f pearMix_hard


