VERSION = c36_jul22

.PHONY: clean clean_$(VERSION) mostlyclean mostlyclean_$(VERSION) current all diff diff_$(VERSION) new

all : current

current : current_$(VERSION)

clean : clean_$(VERSION)

mostlyclean : mostlyclean_$(VERSION)

new : toppar_$(VERSION)

diff : diff_$(VERSION)

current_$(VERSION) : 
	cp -rPp toppar_$(VERSION)/* .
	touch current_$(VERSION)

clean_$(VERSION) :
	rm -rf toppar_$(VERSION).tgz toppar_$(VERSION)_source toppar_$(VERSION)

mostlyclean_$(VERSION) :
	rm -rf toppar_$(VERSION)_source toppar_$(VERSION)

diff_$(VERSION) : toppar_$(VERSION)/*
	$(foreach d,$^,diff -rq $(d) $(notdir $(d));)

toppar_$(VERSION).tgz :
	wget http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/$@ -O $@

toppar_$(VERSION)_source : toppar_$(VERSION).tgz
	mkdir $@ && tar -xf $< -C $@ --strip-components=1

toppar_$(VERSION) : toppar_$(VERSION)_source
	mkdir $@ && cp -rPp $</stream $@/str
	cp -rPp $</toppar_water_ions.str $@/str
	mkdir $@/top &&	cp -rPp $</*.rtf $@/top
	mkdir $@/par &&	cp -rPp $</*.prm $@/par