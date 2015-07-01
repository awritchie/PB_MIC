NAME := PB_MIC
EXTRA_SRC := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(EXTRA_SRC))

#GRVY_INC := $(TACC_GRVY_INC)
#GRVY_LIB := $(TACC_GRVY_LIB)
GRVY_INC := /usr/local/include
GRVY_LIB := /usr/local/lib

CC := g++
INCLUDE := -I. -I$(GRVY_INC) 
LDFLAGS := -L$(GRVY_LIB) 
LDLIBS := -lgrvy -llapack -lblas 
CPPFLAGS := -Wall -O3 $(INCLUDE) -fopenmp

$(NAME): $(OBJS)
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) -o $@ $^
%.o: %.cpp
	$(CC) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS) -c $<
clean:
	rm -f $(NAME) $(OBJS)

.PHONY: clean

