SHELL=/bin/sh
#
SRCS = main.c init.c calc.c cell.c mc.c inout.c function.c ran.c 

OBJS = main.o init.o calc.o cell.o mc.c inout.o function.o ran.o


HFILES =  dpd.h
CF     =  gcc
FFLAGS = 
LIBS   = -L. -lm 
OUT    =  dpd.run 


opt1:	$(SRCS:.c=.o)
	$(CF) $(FFLAGS) -o $(OUT) $(SRCS:.c=.o) $(LIBS)

opt2:	$(SRCS)
	$(CF) $(FFLAGS) -O2 -o $(OUT) $(SRCS) $(LIBS)

opt3:   $(SRCS)
	$(CF) $(FFLAGS) -O3 -o $(OUT) $(SRCS) $(LIBS)

prof:   $(SRCS)
	$(CF) $(FFLAGS) -pg -o $(OUT) $(SRCS) $(LIBS)

$(SRCS):	$(HFILES)
$(OBJS):	$(HFILES)	







