# Example call: (will write to htmd/lib/Linux/pro/)
# make CC=g++ PLATFORM=Linux TYPE=pro

CC := $(CC)
FC := $(FC)
PLATFORM := $(PLATFORM)
TYPE := $(TYPE)

_SRC:= TMalign

_OBJ := $(patsubst %,obj/$(PLATFORM)/%.o, $(_SRC))
_SRC := $(patsubst %,src/$(PLATFORM)/%.cpp, $(_SRC))

LDFLAGS :=
INC_FLAGS :=

.PHONY: clean

all: obj/$(PLATFORM) ../../htmd/lib/$(PLATFORM)/tmalign.so

../../htmd/lib/$(PLATFORM)/tmalign.so: $(_OBJ)
	$(CC) -fPIC  --shared -o $@  $(_OBJ) -lm 

obj/$(PLATFORM)/%.o: %.cpp
	$(CC) -o $@ -c $< -g -fPIC -O3 -DPLATFORM_$(PLATFORM)

../../htmd/lib/$(PLATFORM):
	mkdir -p ../../htmd/lib/$(PLATFORM)

obj/$(PLATFORM):
	mkdir -p obj/$(PLATFORM)

clean:
	rm -rf obj/* obj

exe: TMalign.cpp basic_fun.h  Kabsch.h  NW.h TMalign.h
	$(CC) -static -O3 -ffast-math -lm -o TMalign TMalign.cpp


