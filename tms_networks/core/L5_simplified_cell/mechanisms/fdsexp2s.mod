COMMENT
Implementation of the model of short-term facilitation and depression described in
  Varela, J.A., Sen, K., Gibson, J., Fost, J., Abbott, L.R., and Nelson, S.B.
  A quantitative description of short-term plasticity at excitatory synapses
  in layer 2/3 of rat primary visual cortex
  Journal of Neuroscience 17:7926-7940, 1997
This is a modification of Exp2Syn that can receive multiple streams of
synaptic input via NetCon objects.  Each stream keeps track of its own
weight and activation history.

The printf() statements are for testing purposes only.


The synaptic mechanism itself uses a two state kinetic scheme described by
rise time tau1 and decay time constant tau2.
The normalized peak condunductance is 1.
Decay time MUST be greater than rise time.

The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is
 A = a*exp(-t/tau1) and
 G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))
	where tau1 < tau2

If tau2-tau1 -> 0 then we have a alphasynapse.
and if tau1 -> 0 then we have just single exponential decay.

The factor is evaluated in the
initial block such that an event of weight 1 generates a
peak conductance of 1.

Because the solution is a sum of exponentials, the
coupled equations can be solved as a pair of independent equations
by the more efficient cnexp method.

ENDCOMMENT

NEURON {
	POINT_PROCESS FDSExp2Syn
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i

	RANGE g
	GLOBAL total
        RANGE f, tau_F, d1, tau_D1, d2, tau_D2
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	tau1 = 0.1 (ms) < 1e-9, 1e9 >
	tau2 = 10 (ms) < 1e-9, 1e9 >
	e = 0	(mV)
        : these values are from Fig.3 in Varela et al. 1997
	: the (1) is needed for the range limits to be effective
        f = 0.917 (1) < 0, 1e9 >    : facilitation
        tau_F = 94 (ms) < 1e-9, 1e9 >
        d1 = 0.416 (1) < 0, 1 >     : fast depression
        tau_D1 = 380 (ms) < 1e-9, 1e9 >
        d2 = 0.975 (1) < 0, 1 >     : slow depression
        tau_D2 = 9200 (ms) < 1e-9, 1e9 >
}

ASSIGNED {
	v (mV)
	i (nA)
	g (umho)
	factor
	total (umho)
}

STATE {
	A (umho)
	B (umho)
}

INITIAL {
	LOCAL tp
	total = 0
	if (tau1/tau2 > 0.9999) {
		tau1 = 0.9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = B - A
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

NET_RECEIVE(weight (umho), F, D1, D2, tsyn (ms)) {
INITIAL {
: these are in NET_RECEIVE to be per-stream
        F = 1
        D1 = 1
        D2 = 1
        tsyn = t
: this header will appear once per stream
: printf("t\t t-tsyn\t F\t D1\t D2\t amp\t newF\t newD1\t newD2\n")
}

        F = 1 + (F-1)*exp(-(t - tsyn)/tau_F)
        D1 = 1 - (1-D1)*exp(-(t - tsyn)/tau_D1)
        D2 = 1 - (1-D2)*exp(-(t - tsyn)/tau_D2)
: printf("%g\t%g\t%g\t%g\t%g\t%g", t, t-tsyn, F, D1, D2, weight*F*D1*D2)
        tsyn = t

	state_discontinuity(A, A + weight*factor*F*D1*D2)
	state_discontinuity(B, B + weight*factor*F*D1*D2)
	total = total+weight*F*D1*D2

        F = F + f
        D1 = D1 * d1
        D2 = D2 * d2
: printf("\t%g\t%g\t%g\n", F, D1, D2)
}
