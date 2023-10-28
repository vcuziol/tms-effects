/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__ASyN_STD
#define _nrn_initial _nrn_initial__ASyN_STD
#define nrn_cur _nrn_cur__ASyN_STD
#define _nrn_current _nrn_current__ASyN_STD
#define nrn_jacob _nrn_jacob__ASyN_STD
#define nrn_state _nrn_state__ASyN_STD
#define _net_receive _net_receive__ASyN_STD 
#define state state__ASyN_STD 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define tau1 _p[0]
#define tau1_columnindex 0
#define tau2 _p[1]
#define tau2_columnindex 1
#define e _p[2]
#define e_columnindex 2
#define f _p[3]
#define f_columnindex 3
#define tau_F _p[4]
#define tau_F_columnindex 4
#define d1 _p[5]
#define d1_columnindex 5
#define tau_D1 _p[6]
#define tau_D1_columnindex 6
#define d2 _p[7]
#define d2_columnindex 7
#define tau_D2 _p[8]
#define tau_D2_columnindex 8
#define gpeak _p[9]
#define gpeak_columnindex 9
#define i _p[10]
#define i_columnindex 10
#define g _p[11]
#define g_columnindex 11
#define A _p[12]
#define A_columnindex 12
#define B _p[13]
#define B_columnindex 13
#define factor _p[14]
#define factor_columnindex 14
#define DA _p[15]
#define DA_columnindex 15
#define DB _p[16]
#define DB_columnindex 16
#define _g _p[17]
#define _g_columnindex 17
#define _tsav _p[18]
#define _tsav_columnindex 18
#define _nd_area  *_ppvar[0]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(Object* _ho) { void* create_point_process(int, Object*);
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(void*);
 static double _hoc_loc_pnt(void* _vptr) {double loc_point_process(int, void*);
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(void* _vptr) {double has_loc_point(void*);
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(void* _vptr) {
 double get_loc_point_process(void*); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 0, 0
};
 /* declare global and static user variables */
#define total total_ASyN_STD
 double total = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "d2", 0, 1,
 "d1", 0, 1,
 "f", 0, 1e+09,
 "tau_D2", 1e-09, 1e+09,
 "tau_D1", 1e-09, 1e+09,
 "tau_F", 1e-09, 1e+09,
 "tau2", 1e-09, 1e+09,
 "tau1", 1e-09, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "total_ASyN_STD", "umho",
 "tau1", "ms",
 "tau2", "ms",
 "e", "mV",
 "f", "1",
 "tau_F", "ms",
 "d1", "1",
 "tau_D1", "ms",
 "d2", "1",
 "tau_D2", "ms",
 "gpeak", "umho",
 "A", "umho",
 "B", "umho",
 "i", "nA",
 "g", "umho",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double delta_t = 0.01;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "total_ASyN_STD", &total_ASyN_STD,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"ASyN_STD",
 "tau1",
 "tau2",
 "e",
 "f",
 "tau_F",
 "d1",
 "tau_D1",
 "d2",
 "tau_D2",
 "gpeak",
 0,
 "i",
 "g",
 0,
 "A",
 "B",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 19, _prop);
 	/*initialize range parameters*/
 	tau1 = 0.1;
 	tau2 = 10;
 	e = 0;
 	f = 0.917;
 	tau_F = 50;
 	d1 = 0.5;
 	tau_D1 = 200;
 	d2 = 0.975;
 	tau_D2 = 9200;
 	gpeak = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 19;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _net_receive(Point_process*, double*, double);
 static void _net_init(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _alpha_STD_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 19, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_init[_mechtype] = _net_init;
 pnt_receive_size[_mechtype] = 5;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ASyN_STD alpha_STD.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[2], _dlist1[2];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   DA = - A / tau1 ;
   DB = - B / tau2 ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 DA = DA  / (1. - dt*( ( - 1.0 ) / tau1 )) ;
 DB = DB  / (1. - dt*( ( - 1.0 ) / tau2 )) ;
  return 0;
}
 /*END CVODE*/
 static int state () {_reset=0;
 {
    A = A + (1. - exp(dt*(( - 1.0 ) / tau1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - A) ;
    B = B + (1. - exp(dt*(( - 1.0 ) / tau2)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - B) ;
   }
  return 0;
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{    _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
   _args[2] = 1.0 - ( 1.0 - _args[2] ) * exp ( - ( t - _args[4] ) / tau_D1 ) ;
   _args[4] = t ;
       if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = A;
    double __primary = (A + _args[0] * factor * _args[1] * _args[2] * _args[3] ) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau1 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau1 ) - __primary );
    A += __primary;
  } else {
 A = A + _args[0] * factor * _args[1] * _args[2] * _args[3]  ;
     }
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B;
    double __primary = (B + _args[0] * factor * _args[1] * _args[2] * _args[3] ) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau2 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau2 ) - __primary );
    B += __primary;
  } else {
 B = B + _args[0] * factor * _args[1] * _args[2] * _args[3]  ;
     }
 total = total + _args[0] * _args[1] * _args[2] * _args[3] ;
   _args[2] = _args[2] * d1 ;
   } }
 
static void _net_init(Point_process* _pnt, double* _args, double _lflag) {
       _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
 _args[1] = 1.0 ;
   _args[2] = 1.0 ;
   _args[3] = 1.0 ;
   _args[4] = t ;
   }
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  A = A0;
  B = B0;
 {
   double _ltp ;
 total = 0.0 ;
   if ( tau1 / tau2 > 0.9999 ) {
     tau1 = 0.9999 * tau2 ;
     }
   A = 0.0 ;
   B = 0.0 ;
   _ltp = ( tau1 * tau2 ) / ( tau2 - tau1 ) * log ( tau2 / tau1 ) ;
   factor = - exp ( - _ltp / tau1 ) + exp ( - _ltp / tau2 ) ;
   factor = 1.0 / factor ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   g = gpeak * ( B - A ) ;
   i = g * ( v - e ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 { error =  state();
 if(error){fprintf(stderr,"at line 95 in file alpha_STD.mod:\n	SOLVE state METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 }}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = A_columnindex;  _dlist1[0] = DA_columnindex;
 _slist1[1] = B_columnindex;  _dlist1[1] = DB_columnindex;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "alpha_STD.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "Implementation of the model of short-term facilitation and depression described in\n"
  "  Varela, J.A., Sen, K., Gibson, J., Fost, J., Abbott, L.R., and Nelson, S.B.\n"
  "  A quantitative description of short-term plasticity at excitatory synapses\n"
  "  in layer 2/3 of rat primary visual cortex\n"
  "  Journal of Neuroscience 17:7926-7940, 1997\n"
  "This is a modification of Exp2Syn that can receive multiple streams of\n"
  "synaptic input via NetCon objects.  Each stream keeps track of its own\n"
  "weight and activation history.\n"
  "\n"
  "The printf() statements are for testing purposes only.\n"
  "\n"
  "\n"
  "The synaptic mechanism itself uses a two state kinetic scheme described by\n"
  "rise time tau1 and decay time constant tau2.\n"
  "The normalized peak conductance is 1.\n"
  "Decay time MUST be greater than rise time.\n"
  "\n"
  "The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is\n"
  " A = a*exp(-t/tau1) and\n"
  " G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))\n"
  "	where tau1 < tau2\n"
  "\n"
  "If tau2-tau1 -> 0 then we have a alphasynapse.\n"
  "and if tau1 -> 0 then we have just single exponential decay.\n"
  "\n"
  "The factor is evaluated in the\n"
  "initial block such that an event of weight 1 generates a\n"
  "peak conductance of 1.\n"
  "\n"
  "Because the solution is a sum of exponentials, the\n"
  "coupled equations can be solved as a pair of independent equations\n"
  "by the more efficient cnexp method.\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS ASyN_STD\n"
  "	RANGE tau1, tau2, e, i, gpeak\n"
  "	NONSPECIFIC_CURRENT i\n"
  "\n"
  "	RANGE g\n"
  "	GLOBAL total\n"
  "        RANGE f, tau_F, d1, tau_D1, d2, tau_D2\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	tau1 = 0.1 (ms) < 1e-9, 1e9 >\n"
  "	tau2 = 10 (ms) < 1e-9, 1e9 >\n"
  "	e = 0	(mV)\n"
  "        : these values are from Fig.3 in Varela et al. 1997\n"
  "	: the (1) is needed for the range limits to be effective\n"
  "        f = 0.917 (1) < 0, 1e9 >    : facilitation\n"
  "        tau_F = 50 (ms) < 1e-9, 1e9 >\n"
  "        d1 = 0.5 (1) < 0, 1 >     : fast depression\n"
  "        tau_D1 = 200 (ms) < 1e-9, 1e9 >\n"
  "        d2 = 0.975 (1) < 0, 1 >     : slow depression\n"
  "        tau_D2 = 9200 (ms) < 1e-9, 1e9 >\n"
  "		gpeak		(umho)		: maximum conductance\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	i (nA)\n"
  "	g (umho)\n"
  "	factor\n"
  "	total (umho)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	A (umho)\n"
  "	B (umho)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	LOCAL tp\n"
  "	total = 0\n"
  "	if (tau1/tau2 > 0.9999) {\n"
  "		tau1 = 0.9999*tau2\n"
  "	}\n"
  "	A = 0\n"
  "	B = 0\n"
  "	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)\n"
  "	factor = -exp(-tp/tau1) + exp(-tp/tau2)\n"
  "	factor = 1/factor\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "	g = gpeak *(B - A)\n"
  "	i = g*(v - e)\n"
  "}\n"
  "\n"
  "DERIVATIVE state {\n"
  "	A' = -A/tau1\n"
  "	B' = -B/tau2\n"
  "}\n"
  "\n"
  "NET_RECEIVE(weight (umho), F, D1, D2, tsyn (ms)) {\n"
  "INITIAL {\n"
  ": these are in NET_RECEIVE to be per-stream\n"
  "        F = 1\n"
  "        D1 = 1\n"
  "        D2 = 1\n"
  "        tsyn = t\n"
  ": this header will appear once per stream\n"
  ": printf(\"t\\t t-tsyn\\t F\\t D1\\t D2\\t amp\\t newF\\t newD1\\t newD2\\n\")\n"
  "}\n"
  "\n"
  "        :F = 1 + (F-1)*exp(-(t - tsyn)/tau_F)\n"
  "        D1 = 1 - (1-D1)*exp(-(t - tsyn)/tau_D1)\n"
  "        :D2 = 1 - (1-D2)*exp(-(t - tsyn)/tau_D2)\n"
  ": printf(\"%g\\t%g\\t%g\\t%g\\t%g\\t%g\", t, t-tsyn, F, D1, D2, weight*F*D1*D2)\n"
  "        tsyn = t\n"
  "\n"
  "	state_discontinuity(A, A + weight*factor*F*D1*D2)\n"
  "	state_discontinuity(B, B + weight*factor*F*D1*D2)\n"
  "	total = total+weight*F*D1*D2\n"
  "\n"
  "        :F = F + f\n"
  "        D1 = D1 * d1\n"
  "        :D2 = D2 * d2\n"
  ": printf(\"\\t%g\\t%g\\t%g\\n\", F, D1, D2)\n"
  "}\n"
  ;
#endif
