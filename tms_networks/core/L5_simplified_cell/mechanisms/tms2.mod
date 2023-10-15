COMMENT
Monophasic current induced by a TMS stimulation. Positive values of i depolarize the cell.
ENDCOMMENT

NEURON {
    POINT_PROCESS TMSpulse_bi
    RANGE onset,decay,imax
    ELECTRODE_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (uS) = (microsiemens)
}

PARAMETER {
    onset = 10 (ms)
    decay = 0.08   (ms)
    imax = 40 	(nA)
}

ASSIGNED {
    i (nA)
}

INITIAL {
	i = 40 (nA)
}
BREAKPOINT {
	if (imax) {
	at_time(onset)
	}

	if ((t - onset) < 0 || (t - onset) > 20) {
	   i = 0
	 } else {
	 i = imax*sin (30*(t - onset)) * exp ( -(t - onset)/decay)
	 }


}

FUNCTION alpha(x) {
	if (x < 0 || x > 20) {
		alpha = 0
	}else{
		alpha = exp(1 - x) - x * exp(1 - x)
	}
}
