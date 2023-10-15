TITLE HH channel
: Mel-modified Hodgkin - Huxley conductances (after Ojvind et al.)

: modified by C. Rusu in 2012 (removed variable s)

VERBATIM
static const char rcsid[]="$Id: hh3.mod,v 1.1 1996/05/19 19:26:28 karchie Exp $";
ENDVERBATIM

NEURON {
	SUFFIX hh3
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gkbar, gl, el,gkbar2,vshift
	GLOBAL taus,taun,taum,tauh,tausb,taun2
	GLOBAL tausv,tausd,mN,nN,sN
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	celsius = 37	(degC)
	dt (ms)
	gnabar=.20 (mho/cm2)
	gkbar=.12 (mho/cm2)
	gkbar2=.12 (mho/cm2)
	gl=.0001 (mho/cm2)
	ena = 40 (mV)
	ek = -80 (mV)
	el = -70.0 (mV)	: steady state at v = -65 mV
	taum=0.05
	tauh=0.1 : originally 0.5
	taus=50
	tausv=30
	tausd=1
	taun=1
	taun2	=10
	mN=3
	nN=3
	sN=0            : originally 1
	tausb=0.5
	vshift=0
}
STATE {
	m h n s n2
}
ASSIGNED {
	ina (mA/cm2)
	ik (mA/cm2)
	il (mA/cm2)

}

BREAKPOINT {
	SOLVE states

	ina = gnabar*h*s^sN*(v - ena)*m^mN
	ik = gkbar*(v - ek)*n^nN+gkbar2*(v - ek)*n2^nN


	il = gl*(v - el)
}

PROCEDURE states() {	: exact when v held constant
	LOCAL sigmas
	sigmas=1/(1+exp((v+tausv+vshift)/tausd))
	m = m + (1 - exp(-dt/taum))*(1 / (1 + exp((v + 40+vshift)/(-3)))  - m)
	h = h + (1 - exp(-dt/tauh))*(1 / (1 + exp((v + 45+vshift)/3))  - h)
	s = s + (1 - exp(-dt/(taus*sigmas+tausb)))*(1 / (1 + exp((v + 44+vshift)/3))  - s)
	n = n + (1 - exp(-dt/taun))*(1 / (1 + exp((v + 40+vshift)/(-3)))  - n)
	n2 = n2 + (1 - exp(-dt/taun2))*(1 / (1 + exp((v + 40+vshift)/(-3)))  - n2)
	VERBATIM
	return 0;
	ENDVERBATIM
}

