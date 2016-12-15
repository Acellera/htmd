CC := $(CC)
CPURE := $(CPURE)
FC := $(FC)
PLATFORM := $(PLATFORM)
TYPE := $(TYPE)

SRC:= occupancy_ext

_OBJ := $(patsubst %,obj/$(PLATFORM)/%.o, $(SRC))
_SRC := $(patsubst %,src/$(PLATFORM)/%.c, $(SRC))

LDFLAGS := 
INC_FLAGS := 

.PHONY: clean 

all: obj/$(PLATFORM) ../../htmd/lib/$(PLATFORM)/occupancy_ext.so

../../htmd/lib/$(PLATFORM)/occupancy_ext.so: $(_OBJ)
	$(CPURE) -fPIC  --shared -o $@  $(_OBJ)

obj/$(PLATFORM)/%.o: %.c
	$(CPURE) -o $@ -c $< -g -fPIC -ffast-math -O3 -DDIST_

../../htmd/lib/$(PLATFORM):
	mkdir -p ../../htmd/lib/$(PLATFORM)

obj/$(PLATFORM):
	mkdir -p obj/$(PLATFORM)

clean:
	rm -rf obj/* obj


