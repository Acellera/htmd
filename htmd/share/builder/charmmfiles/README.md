# Notes

How to use:
Update the first line in the Makefile to the latest CHARMM36 version. Then:

```
make new
make diff
```

When doing a `make diff`, the following differences are expected:

```
Only in str/lipid: toppar_all36_lipid_cholesterol_model_1.str
Only in top: top_all22star_prot.rtf
Only in par: par_all22star_prot.prm
```

The `toppar_all36_lipid_cholesterol_model_1.str` is a processed file created by us to use the right cholesterol.
The `top_all22star_prot.rtf` and `par_all22star_prot.prm` are files added by us that are not present in the standard CHARMM force-field.

Manually overwrite the old parameters with the new ones.
