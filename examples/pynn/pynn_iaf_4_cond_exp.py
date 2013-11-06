import pyNN.nest as sim
sim.nest.Install('coronetmodule')
sim.setup()
sim.nest.SetKernelStatus({'dict_miss_is_error': False})

iaf_4_cond_exp = sim.native_cell_type('iaf_4_cond_exp')
cell_params=\
{'C_m': 250.0,
 'E_L': -70.0,
 'E_syn_1': 0.0,
 'E_syn_2': 0.0,
 'E_syn_3': -85.0,
 'E_syn_4': 0.0,
 'I_e': 0.0,
 'V_reset': -60.0,
 'V_th': -55.0,
 'g_L': 16.6667,
 't_ref': 2.0,
 'tau_syn_1': 5.0,
 'tau_syn_2': 5.0,
 'tau_syn_3': 25.0,
 'tau_syn_4': 0.2}
p1 = sim.Population(10, iaf_4_cond_exp, cell_params)
s1 = sim.Population(10, sim.SpikeSourcePoisson, {'rate':20})
s2 = sim.Population(10, sim.SpikeSourcePoisson, {'rate':20})
s3 = sim.Population(10, sim.SpikeSourcePoisson, {'rate':20})
sim.Projection(s1, p1, sim.OneToOneConnector(weights=0.02), target="SYN_1")
sim.Projection(s2, p1, sim.OneToOneConnector(weights=0.01), target="SYN_2")
sim.Projection(s3, p1, sim.OneToOneConnector(weights=0.01), target="SYN_3")

print iaf_4_cond_exp.synapse_types

p1.initialize('V_m', -76.) # this works as expected

p1.record()
p1._record('V_m') # ugly
p1._record('g_syn_1') # ugly
p1._record('g_syn_2') # ugly
p1._record('g_syn_3') # ugly
sim.run(1000)
id, t, v = p1.recorders['V_m'].get().T # ugly
id, t, g1 = p1.recorders['g_syn_1'].get().T # ugly
id, t, g2 = p1.recorders['g_syn_2'].get().T # ugly
id, t, g3 = p1.recorders['g_syn_3'].get().T # ugly

import pylab as pl
import numpy as np
pl.figure()
sl = p1.getSpikes()
pl.plot(sl[:,1],sl[:,0],'.')
pl.figure()
id_is_0 = np.where(id==0)
pl.plot(t[id_is_0],v[id_is_0])
pl.show()
pl.plot(t[id_is_0],g1[id_is_0])
pl.plot(t[id_is_0],g2[id_is_0])
pl.plot(t[id_is_0],g3[id_is_0])
pl.show()
