VERSION = c36_jul18

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
	mkdir $@/top &&	cp -rPp $</*.rtf $@/top
	mkdir $@/par &&	cp -rPp $</*.prm $@/par
# hacky split the water file
	awk '/read rtf card \@app/,/END/{if (!/read rtf card \@app/){print}}' $</toppar_water_ions.str > $@/top/top_water_ions.rtf
	awk '/read para card flex \@app/,NR==325{if (!/read para card flex \@app/){print}}' $</toppar_water_ions.str > $@/par/par_water_ions.prm
