CXX=gcc
LIBS=-lgslcblas -lgsl -lblas -L/usr/lib -lm  -L/usr/lib
CFLAGS=-fPIC
CFLAGS_DEBUG=-g -pg -Wall -W

LIBFILE=sct_smart_boxes.so
OBJS=sct_smart_boxes.o

ifneq (,$(findstring debug,$(MAKECMDGOALS)))
  override CFLAGS += $(CFLAGS_DEBUG)
endif

all: $(LIBFILE)
debug: all

$(LIBFILE): $(OBJS)
	$(CXX) $< $(CFLAGS) $(LIBS) -shared -O3 -o $(LIBFILE)

%.o : %.c
	$(CXX) $< $(CFLAGS) -c $(LIBS) -o $@

clean:
	rm -f $(LIBFILE) $(OBJS)
