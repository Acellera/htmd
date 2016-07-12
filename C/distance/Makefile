CC := $(CC)
FC := $(FC)
PLATFORM := $(PLATFORM)
TYPE := $(TYPE)

MINDIST_SRC:= mindist_ext

_MINDIST_OBJ := $(patsubst %,obj/$(TYPE)/$(PLATFORM)/%.o, $(MINDIST_SRC))
_MINDIST_SRC := $(patsubst %,src/$(PLATFORM)/%.c, $(MINDIST_SRC))

LDFLAGS := 
INC_FLAGS := 

.PHONY: clean 

dummy:
	./build.sh

all: obj/$(TYPE)/$(PLATFORM) ../../htmd/lib/$(TYPE)/$(PLATFORM)/mindist_ext.so

../../htmd/lib/$(TYPE)/$(PLATFORM)/mindist_ext.so: $(_MINDIST_OBJ)
	$(CC) -fPIC  --shared -o $@  $(_MINDIST_OBJ)

obj/$(TYPE)/$(PLATFORM)/%.o: %.c
	$(CC) -o $@ -c $< -g -fPIC -ffast-math -O3 -DDIST_$(TYPE)

../../htmd/lib/$(TYPE)/$(PLATFORM):
	mkdir -p ../../htmd/lib/$(TYPE)/$(PLATFORM)

obj/$(TYPE)/$(PLATFORM):
	mkdir -p obj/$(TYPE)/$(PLATFORM)

clean:
	rm -rf obj/* obj


