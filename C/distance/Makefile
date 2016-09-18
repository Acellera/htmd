CPURE := $(CPURE)
FC := $(FC)
PLATFORM := $(PLATFORM)
TYPE := $(TYPE)

MINDIST_SRC:= mindist_ext

_MINDIST_OBJ := $(patsubst %,obj/$(PLATFORM)/%.o, $(MINDIST_SRC))
_MINDIST_SRC := $(patsubst %,src/$(PLATFORM)/%.c, $(MINDIST_SRC))

LDFLAGS := 
INC_FLAGS := 

.PHONY: clean 

all: obj/$(PLATFORM) ../../htmd/lib/$(PLATFORM)/mindist_ext.so

../../htmd/lib/$(PLATFORM)/mindist_ext.so: $(_MINDIST_OBJ)
	$(CPURE) -fPIC  --shared -o $@  $(_MINDIST_OBJ)

obj/$(PLATFORM)/%.o: %.c
	$(CPURE) -o $@ -c $< -g -fPIC -ffast-math -O3 -DDIST_

../../htmd/lib/$(PLATFORM):
	mkdir -p ../../htmd/lib/$(PLATFORM)

obj/$(PLATFORM):
	mkdir -p obj/$(PLATFORM)

clean:
	rm -rf obj/* obj


