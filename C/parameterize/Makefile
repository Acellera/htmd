CC := $(CC) -static
FC := $(FC) -static

TIMEBOMB := -DTIMEBOMB=1 -DCDATE=$(shell date +%s)
CFLAGS := -g -Wno-format-zero-length -Wno-unused-but-set-variable -Wno-unused-variable -Wno-write-strings -Wall #-fsanitize=address 
LDFLAGS:= -Lnlopt/lib -lnlopt

NLOPT  := nlopt/lib/libnlopt.a

CPROGS := 1d_fitting 1d_org 1d_rotamer_fitting org_rotamer_E acceptor  \
         add_tip3 check_b0_theta0 check_lj clustering_phi gen_soft_list donor \
         equiv_atom fitcharge fitcharge_again gen_esp gen_xpsf mm_pes \
         mm_pes_large qm_1d_scan_large_para qm_1d_scan_para \
         qm_1d_scan_single_para qm_rotamer_scan qm_rotamer_scan_large \
         update_tor_para update_xpsf  don_acc_separate  \
				 1d_fitting_asymmetric

#drude_1d_fitting \
#         drude_gen_esp drude_gen_xpsf drude_modify_lj drude_ff drude_scale_alpha \
#         drude_scale_vdw_r6 drude_update_torsion_para

DEST=../../htmd/parameterize/share/bin
CBIN := $(patsubst %,$(DEST)/%, $(CPROGS))
FBIN := $(DEST)/cgrid # $(DEST)/drude_cgrid

.PHONY: clean
.PHONY: scripts

all: $(DEST) lib lib/libff.a lib/libdrudeff.a $(NLOPT) $(CBIN) $(FBIN) scripts


scripts:
	cp scripts/* $(DEST)/

$(DEST)/drude_cgrid:
	$(FC) -static -o $@ src/drude_ConnollyGrid.f src/drude_Surface.f

$(DEST)/cgrid:
	$(FC) -static -o $@  src/ConnollyGrid.f src/Surface.f

$(DEST)/drude_fitcharge: src/drude-fitcharge.cpp
	$(CC) -static -o $@ $< $(CFLAGS) $(LDFLAGS) -Llib -ldrude-ff

$(DEST)/drude_fitcharge_again: src/drude-fitcharge_again.cpp
	$(CC) -static -o $@ $< $(CFLAGS) $(LDFLAGS) -Llib -ldrude-ff

$(DEST)/%: src/%.cpp lib/libff.a
	$(CC) -static -o $@  $< $(CFLAGS) $(LDFLAGS) -Llib -lff -lm



lib/libff.a:  src/ff.cpp src/ff.h src/timebomb.cpp
	$(CC) $(CFLAGS) -o $@  -c src/ff.cpp $(TIMEBOMB)

lib/libdrudeff.a: 
	$(CC) $(CFLAGS) -o $@ -c src/ff.cpp $(TIMEBOMB)


$(NLOPT):
	nlopt/build-nlopt.sh

clean:
	rm -rf lib/* $(DEST)/* nlopt/build nlopt/lib nlopt/include  

lib:
	mkdir -p lib

$(DEST):
	mkdir -p $@
