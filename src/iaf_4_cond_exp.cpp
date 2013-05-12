/*
 *  iaf_4_cond_exp.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "iaf_4_cond_exp.h"

#ifdef HAVE_GSL

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include <limits>

#include "universal_data_logger_impl.h"
#include "event.h"

#include <iomanip>
#include <iostream>
#include <cstdio>

using namespace nest;

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<mynest::iaf_4_cond_exp> mynest::iaf_4_cond_exp::recordablesMap_;

namespace nest  // template specialization must be placed in namespace
{
  // Override the create() method with one call to RecordablesMap::insert_() 
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<mynest::iaf_4_cond_exp>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(names::V_m, 
	    &mynest::iaf_4_cond_exp::get_y_elem_<mynest::iaf_4_cond_exp::State_::V_M>);
    insert_("g_syn_1", 
	    &mynest::iaf_4_cond_exp::get_y_elem_<mynest::iaf_4_cond_exp::State_::G_SYN_1>);
    insert_("g_syn_2", 
	    &mynest::iaf_4_cond_exp::get_y_elem_<mynest::iaf_4_cond_exp::State_::G_SYN_2>);
    insert_("g_syn_3", 
	    &mynest::iaf_4_cond_exp::get_y_elem_<mynest::iaf_4_cond_exp::State_::G_SYN_3>);
    insert_("g_syn_4", 
	    &mynest::iaf_4_cond_exp::get_y_elem_<mynest::iaf_4_cond_exp::State_::G_SYN_4>);
  }
}

extern "C"
inline int mynest::iaf_4_cond_exp_dynamics(double, const double y[], double f[], void* pnode)
{ 
  // a shorthand
  typedef mynest::iaf_4_cond_exp::State_ S;

  // get access to node so we can almost work as in a member function
  assert(pnode);
  const mynest::iaf_4_cond_exp& node =  *(reinterpret_cast<mynest::iaf_4_cond_exp*>(pnode));

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[]. 

  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...
  // TODO: make code nicer through looping
  const double I_syn_1 = y[S::G_SYN_1] * (y[S::V_M] - node.P_.E_syn_1); 
  const double I_syn_2 = y[S::G_SYN_2] * (y[S::V_M] - node.P_.E_syn_2); 
  const double I_syn_3 = y[S::G_SYN_3] * (y[S::V_M] - node.P_.E_syn_3); 
  const double I_syn_4 = y[S::G_SYN_4] * (y[S::V_M] - node.P_.E_syn_4); 
  const double I_L       = node.P_.g_L * (y[S::V_M] - node.P_.E_L );

  //V dot
  f[0]= ( - I_L + node.B_.I_stim_ + node.P_.I_e - I_syn_1 - I_syn_2 - I_syn_3 - I_syn_4) / node.P_.C_m;

  f[S::G_SYN_1] = -y[S::G_SYN_1] / node.P_.tau_syn_1;
  f[S::G_SYN_2] = -y[S::G_SYN_2] / node.P_.tau_syn_2;
  f[S::G_SYN_3] = -y[S::G_SYN_3] / node.P_.tau_syn_3;
  f[S::G_SYN_4] = -y[S::G_SYN_4] / node.P_.tau_syn_4;

  return GSL_SUCCESS;
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
mynest::iaf_4_cond_exp::Parameters_::Parameters_()
  : V_th_      (-55.0    ),  // mV
    V_reset_   (-60.0    ),  // mV
    t_ref_     (  2.0    ),  // ms
    g_L        ( 16.6667 ),  // nS
    C_m        (250.0    ),  // pF
    E_syn_1    (  0.0    ),  // mV
    E_syn_2    (  0.0    ),  // mV
    E_syn_3    (  0.0    ),  // mV
    E_syn_4    (  0.0    ),  // mV
    E_L        (-70.0    ),  // mV
    tau_syn_1  (  0.2    ),  // ms
    tau_syn_2  (  0.2    ),  // ms
    tau_syn_3  (  0.2    ),  // ms
    tau_syn_4  (  0.2    ),  // ms
    I_e        (  0.0    )   // pA
{
}

mynest::iaf_4_cond_exp::State_::State_(const Parameters_& p)
  : r_(0)
{
  y_[V_M] = p.E_L;
  y_[G_SYN_1] = 0.;
  y_[G_SYN_2] = 0.;
  y_[G_SYN_3] = 0.;
  y_[G_SYN_4] = 0.;
}

mynest::iaf_4_cond_exp::State_::State_(const State_& s)
  : r_(s.r_)
{
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = s.y_[i];
}

mynest::iaf_4_cond_exp::State_& mynest::iaf_4_cond_exp::State_::operator=(const State_& s)
{
  assert(this != &s);  // would be bad logical error in program
  
  for ( size_t i = 0 ; i < STATE_VEC_SIZE ; ++i )
    y_[i] = s.y_[i];
  r_ = s.r_;
  return *this;
}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void mynest::iaf_4_cond_exp::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d,names::V_th,         V_th_);
  def<double>(d,names::V_reset,      V_reset_);
  def<double>(d,names::t_ref,        t_ref_);
  def<double>(d,names::g_L,          g_L);
  def<double>(d,names::E_L,          E_L); 
  def<double>(d,"E_syn_1",           E_syn_1);
  def<double>(d,"E_syn_2",           E_syn_2);
  def<double>(d,"E_syn_3",           E_syn_3);
  def<double>(d,"E_syn_4",           E_syn_4);
  def<double>(d,names::C_m,          C_m);
  def<double>(d,"tau_syn_1",   tau_syn_1);
  def<double>(d,"tau_syn_2",   tau_syn_2);
  def<double>(d,"tau_syn_3",   tau_syn_3);
  def<double>(d,"tau_syn_4",   tau_syn_4);
  def<double>(d,names::I_e,          I_e);
}

void mynest::iaf_4_cond_exp::Parameters_::set(const DictionaryDatum& d)
{
  // allow setting the membrane potential
  updateValue<double>(d,names::V_th,    V_th_);
  updateValue<double>(d,names::V_reset, V_reset_);
  updateValue<double>(d,names::t_ref,   t_ref_);
  updateValue<double>(d,names::E_L,     E_L);
  
  updateValue<double>(d,"E_syn_1",      E_syn_1);
  updateValue<double>(d,"E_syn_2",      E_syn_2);
  updateValue<double>(d,"E_syn_3",      E_syn_3);
  updateValue<double>(d,"E_syn_4",      E_syn_4);
  
  updateValue<double>(d,names::C_m,     C_m);
  updateValue<double>(d,names::g_L,     g_L);

  updateValue<double>(d,"tau_syn_1",    tau_syn_1);
  updateValue<double>(d,"tau_syn_2",    tau_syn_2);
  updateValue<double>(d,"tau_syn_3",    tau_syn_3);
  updateValue<double>(d,"tau_syn_4",    tau_syn_4);

  updateValue<double>(d,names::I_e,     I_e);

  if ( V_reset_ >= V_th_ )
    throw BadProperty("Reset potential must be smaller than threshold.");
    
  if ( C_m <= 0 )
    throw BadProperty("Capacitance must be strictly positive.");
    
  if ( t_ref_ < 0 )
    throw BadProperty("Refractory time cannot be negative.");
      
  if ( tau_syn_1 <= 0 || tau_syn_2 <= 0 || tau_syn_3 <= 0 || tau_syn_4 <= 0 )
    throw BadProperty("All time constants must be strictly positive.");
}

void mynest::iaf_4_cond_exp::State_::get(DictionaryDatum &d) const
{
  def<double>(d, names::V_m, y_[V_M]); // Membrane potential
  def<double>(d, "g_syn_1", y_[G_SYN_1]);
  def<double>(d, "g_syn_2", y_[G_SYN_2]);
  def<double>(d, "g_syn_3", y_[G_SYN_3]);
  def<double>(d, "g_syn_4", y_[G_SYN_4]);
}

void mynest::iaf_4_cond_exp::State_::set(const DictionaryDatum& d, const Parameters_&)
{
  updateValue<double>(d, names::V_m, y_[V_M]);
  updateValue<double>(d, "g_syn_1", y_[G_SYN_1]);
  updateValue<double>(d, "g_syn_2", y_[G_SYN_2]);
  updateValue<double>(d, "g_syn_3", y_[G_SYN_3]);
  updateValue<double>(d, "g_syn_4", y_[G_SYN_4]);
}

mynest::iaf_4_cond_exp::Buffers_::Buffers_(iaf_4_cond_exp& n)
  : logger_(n),
    spike_inputs_(std::vector<RingBuffer>(SUP_SPIKE_RECEPTOR-1)),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

mynest::iaf_4_cond_exp::Buffers_::Buffers_(const Buffers_&, iaf_4_cond_exp& n)
  : logger_(n),
    spike_inputs_(std::vector<RingBuffer>(SUP_SPIKE_RECEPTOR-1)),
    s_(0),
    c_(0),
    e_(0)
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

mynest::iaf_4_cond_exp::iaf_4_cond_exp()
  : Archiving_Node(), 
    P_(), 
    S_(P_),
    B_(*this)
{
  recordablesMap_.create();
}

mynest::iaf_4_cond_exp::iaf_4_cond_exp(const iaf_4_cond_exp& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
}

mynest::iaf_4_cond_exp::~iaf_4_cond_exp()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
  if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
  if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void mynest::iaf_4_cond_exp::init_state_(const Node& proto)
{
  const iaf_4_cond_exp& pr = downcast<iaf_4_cond_exp>(proto);
  S_ = pr.S_;
}

void mynest::iaf_4_cond_exp::init_buffers_()
{

    // Reset spike buffers.
    for(std::vector<RingBuffer>::iterator it = B_.spike_inputs_.begin();
	it != B_.spike_inputs_.end(); ++it)
      {
	it->clear(); // include resize
      }

  B_.currents_.clear();           // includes resize
  Archiving_Node::clear_history();

  B_.logger_.reset();

  B_.step_ = Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  static const gsl_odeiv_step_type* T1 = gsl_odeiv_step_rkf45;
  
  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc (T1, State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_step_reset(B_.s_);
    
  if ( B_.c_ == 0 )  
    B_.c_ = gsl_odeiv_control_y_new (1e-3, 0.0);
  else
    gsl_odeiv_control_init(B_.c_, 1e-3, 0.0, 1.0, 0.0);
    
  if ( B_.e_ == 0 )  
    B_.e_ = gsl_odeiv_evolve_alloc(State_::STATE_VEC_SIZE);
  else 
    gsl_odeiv_evolve_reset(B_.e_);
  
  B_.sys_.function  = iaf_4_cond_exp_dynamics; 
  B_.sys_.jacobian  = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params    = reinterpret_cast<void*>(this);

  B_.I_stim_ = 0.0;
}

void mynest::iaf_4_cond_exp::calibrate()
{
  B_.logger_.init();  // ensures initialization in case mm connected after Simulate

  V_.RefractoryCounts_ = Time(Time::ms(P_.t_ref_)).get_steps();
  assert(V_.RefractoryCounts_ >= 0);  // since t_ref_ >= 0, this can only fail in error
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void mynest::iaf_4_cond_exp::update(Time const & origin, const long_t from, const long_t to)
{
   
  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);

  for ( long_t lag = from ; lag < to ; ++lag )
  {
    
    double t = 0.0;

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( t < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply(B_.e_, B_.c_, B_.s_, 
			   &B_.sys_,             // system of ODE
			   &t,                   // from t
			    B_.step_,            // to t <= step
			   &B_.IntegrationStep_, // integration step size
			    S_.y_); 	         // neuronal state

      if ( status != GSL_SUCCESS )
        throw GSLSolverFailure(get_name(), status);
    }

    S_.y_[State_::G_SYN_1] += B_.spike_inputs_[Buffers_::SYN_1].get_value(lag);
    S_.y_[State_::G_SYN_2] += B_.spike_inputs_[Buffers_::SYN_2].get_value(lag);
    S_.y_[State_::G_SYN_3] += B_.spike_inputs_[Buffers_::SYN_3].get_value(lag);
    S_.y_[State_::G_SYN_4] += B_.spike_inputs_[Buffers_::SYN_4].get_value(lag);

    // absolute refractory period
    if ( S_.r_ )
    {// neuron is absolute refractory
      --S_.r_; 
      S_.y_[State_::V_M] = P_.V_reset_; 
    }
    else
      // neuron is not absolute refractory
      if ( S_.y_[State_::V_M] >= P_.V_th_ )
	    {
	      S_.r_              = V_.RefractoryCounts_;
	      S_.y_[State_::V_M] = P_.V_reset_;

	      set_spiketime(Time::step(origin.get_steps()+lag+1));
	  
	      SpikeEvent se;
	      network()->send(*this, se, lag);
	    }
    
    // set new input current
    B_.I_stim_ = B_.currents_.get_value(lag);

    // log state data
    B_.logger_.record_data(origin.get_steps() + lag);

  }
}

void mynest::iaf_4_cond_exp::handle(SpikeEvent & e)
{
    assert(e.get_delay() > 0);
    assert(e.get_rport() < static_cast<int_t>(B_.spike_inputs_.size()));

    B_.spike_inputs_[e.get_rport()].
      add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
		e.get_weight() * e.get_multiplicity() );
}

void mynest::iaf_4_cond_exp::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double_t c=e.get_current();
  const double_t w=e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		      w *c);
}

void mynest::iaf_4_cond_exp::handle(DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}

#endif //HAVE_GSL
