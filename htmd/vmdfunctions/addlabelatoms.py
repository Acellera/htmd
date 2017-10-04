

def addLabelAtoms(selection, txtsize=0.5, textcolor='green'):
	return '''label textsize {}
color Labels Atoms {}
set molnum [molinfo top]
set molselected [atomselect $molnum "{}"]
set atomsindexes [$molselected get index]
set i 0
foreach n $atomsindexes {{
	label add Atoms $molnum/$n
	label textformat Atoms $i {{%a}}
	incr i
}} '''.format(txtsize, textcolor, selection,)