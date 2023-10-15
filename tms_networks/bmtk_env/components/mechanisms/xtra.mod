: $Id: xtra.mod,v 1.4 2014/08/18 23:15:25 ted Exp ted $
: 2018/05/20 Modified by Aman Aberra 
: 2019 vcuziol: I removed the 'initial' and 'before breakpoint' sections to be able to manipulate e_extracellular directly.

NEURON {
	SUFFIX xtra
	RANGE es : (es = max amplitude of the potential)		
	RANGE x, y, z, type, order
	GLOBAL stim : (stim = normalized waveform)
	POINTER ex 
}

PARAMETER {	
	es = 0 (mV)
	x = 0 (1) : spatial coords
	y = 0 (1)
	z = 0 (1)		
	type = 0 (1) : numbering system for morphological category of section - unassigned is 0
	order = 0 (1) : order of branch/collateral. 
}

ASSIGNED {
	v (millivolts)
	ex (millivolts)
	stim (unitless) 		
	area (micron2)
}

