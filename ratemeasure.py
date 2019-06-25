#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
import cantera as ct
import timeit
import os
import scipy.optimize as op
import argparse
def runsim ( tmax, temp, pres, ics ):
    dt = tmax/Npoints
    gas.TPX=temp,pres,ics

    r = ct.IdealGasReactor(gas, name='R1')
    sim = ct.ReactorNet([r])

    nr=gas.n_reactions
    for i in range(0,nr):
        r.add_sensitivity_reaction(i)

    sim.rtol = 1.0e-8
    sim.atol = 1.0e-14
    sim.rtol_sensitivity = 1.0e-8
    sim.atol_sensitivity = 1.0e-8

    states = ct.SolutionArray(gas, extra=['t','sens'])
    for i in range(0, int(Npoints)):
        t=i*dt
        sim.advance(t)
        print("%5f"%(t/tmax), end='\r')
        sys.stdout.flush()

        sensitivities=[]
        for i in range(0,nr):
            sensitivities.append(sim.sensitivity('O', i))
        states.append(r.thermo.state, t=t, sens=sensitivities)
    return states

def runsim_nosense ( tmax, temp, pres, ics ):
    dt = tmax/Npoints
    gas.TPX=temp,pres,ics
    r = ct.IdealGasReactor(gas, name='R1')
    sim = ct.ReactorNet([r])

    sim.rtol = 1.0e-8
    sim.atol = 1.0e-14

    states = ct.SolutionArray(gas, extra=['t'])
    for i in range(0, int(Npoints)):
        t=i*dt
        sim.advance(t)

        states.append(r.thermo.state, t=t)

    return states


def residual ( eta, k0, observations, measure_ind, tmaxes, temperatures, pressures, initials, maxes, yields ):
    rstart=timeit.default_timer()

    ret=0;
    reactions = gas.reaction_equations()
    k = k0*10**eta
    if(gas.reaction_type(measure_ind) == 1):
        newrate=ct.ElementaryReaction(gas.reactions()[measure_ind].reactants,gas.reactions()[measure_ind].products)
        newrate.rate=ct.Arrhenius(k, gas.reactions()[measure_ind].rate.temperature_exponent, gas.reactions()[measure_ind].rate.activation_energy)
        gas.modify_reaction(measure_ind, newrate)
    if (gas.reaction_type(measure_ind) == 2):
        newrate=ct.ThreeBodyReaction(reactants=gas.reactions()[measure_ind].reactants,products=gas.reactions()[measure_ind].products)
        newrate.efficiencies = gas.reactions()[measure_ind].efficiencies
        newrate.rate=ct.Arrhenius(k, gas.reactions()[measure_ind].rate.temperature_exponent, gas.reactions()[measure_ind].rate.activation_energy)
        gas.modify_reaction(measure_ind, newrate)
    if (gas.reaction_type(measure_ind) == 4):
        newrate=ct.FalloffReaction(reactants=gas.reactions()[measure_ind].reactants,products=gas.reactions()[measure_ind].products)
        newrate.efficiencies = gas.reactions()[measure_ind].efficiencies
        newrate.falloff = gas.reactions()[measure_ind].falloff
        newrate.low_rate=ct.Arrhenius(k, gas.reactions()[measure_ind].low_rate.temperature_exponent, gas.reactions()[measure_ind].low_rate.temperature_exponent, gas.reactions()[measure_ind].low_rate.activation_energy)
        newrate.high_rate=gas.reactions()[measure_ind].high_rate
        gas.modify_reaction(measure_ind, newrate)

    sys.stdout.flush()
    for i in range(0,nmeasure):
        states = runsim_nosense ( tmaxes[i], temperatures[i], pressures[i], initials[i] )

        for n in maxes:
            lst1=observations[i].X[:,n]
            lst2=states.X[:,n]
            ind1=np.array(np.where(np.sign(np.diff(lst1))==-1))[0,0]

            if np.any(np.sign(np.diff(lst2))==-1):
                ind2=np.array(np.where(np.sign(np.diff(lst2))==-1))[0,0]
            else:
                ind2=len(lst2)-1
            ret+=((1.0*ind2/ind1-1)**2 + (lst2[ind2]/lst1[ind1]-1)**2)/nmeasure

            if(outflag==1):
                print(temperatures[i], tmaxes[i], pressures[i]/ct.one_atm, (1.0*ind2/ind1-1), (lst2[ind2]/lst1[ind1]-1))
                sys.stdout.flush()
        for n in yields:
            lst1=observations[i].X[:,n]
            lst2=states.X[:,n]
            ret+=((lst2[-1]/lst1[-1]-1)**2)/nmeasure

            if(outflag==1):
                print(temperatures[i], tmaxes[i], pressures[i]/ct.one_atm, (1.0*ind2/ind1-1), (lst2[ind2]/lst1[ind1]-1))
                sys.stdout.flush()


    if(outflag==1):
        rstop=timeit.default_timer()
        print ('iter time: %f'%(rstop-rstart))
        print('residual: ', ret)
        print('x: ', eta)
        sys.stdout.flush()

    return np.sqrt(ret)

def remove  ( seed, num, exclude ):
    np.random.seed(seed)
    reactions = gas.reactions()
    removed = []
    rcandidates=range(0,len(reactions))
    rcandidates=[item for item in rcandidates if item not in exclude]
    removed = np.random.choice(rcandidates,size=num,replace=False)
    for ind in removed:
        if(gas.reaction_type(ind) == 1):
            newreac=ct.ElementaryReaction(reactants=reactions[ind].reactants,products=reactions[ind].products)
            gas.modify_reaction(ind, newreac)
        if (gas.reaction_type(ind) == 2):
            newreac=ct.ThreeBodyReaction(reactants=reactions[ind].reactants,products=reactions[ind].products)
            gas.modify_reaction(ind, newreac)
        if (gas.reaction_type(ind) == 4):
            newreac=ct.FalloffReaction(reactants=reactions[ind].reactants,products=reactions[ind].products)
            gas.modify_reaction(ind, newreac)

    return removed

def maxsens (sensitivities):
    maxsens=np.zeros(gas.n_reactions)
    for j in range(0,len(sensitivities)):
        sensitivity=sensitivities[j]
        for i in range(0, sensitivity.shape[1]):
            senslst = (sensitivity.transpose()[i])
            maxsens[i] += (np.linalg.norm(senslst)/len(senslst)**0.5)/len(sensitivities)
    maxind=np.argsort(maxsens)[::-1]
    return maxind, maxsens

parser = argparse.ArgumentParser(description='Measure a rate constant in an incomplete mechanism.')
parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='String for the output files.')
parser.add_argument("--mechanism", type=str, required=False, default='mechanisms/gri30.cti', dest='mechanism', help='Cantera mechanism file')
parser.add_argument("--Npoints", type=float, required=False, default=5e3, dest='Npoints', help='Number of time points to output')
parser.add_argument("--experiments", type=str, required=False, default='experiments/air.dat', dest='experiments', help='File containing a line [tmax temperature pressure initials] for each experimental condition')
parser.add_argument("--measure", type=int, required=False, default=37, dest='measure', help='the index of the reaction whose rate is measured')
parser.add_argument("--remove", type=int, required=False, default=40, dest='remove', help='The number of reactions to randomly remove')
parser.add_argument("--retain", type=int, required=False, default=40, dest='retain')
parser.add_argument("--seed", type=int, required=False, default=1, dest='seed', help='The random seed')
parser.add_argument("--ktol", type=float, required=False, default=1e-6, dest='ktol', help='The tolerance in k relative to k0 in the minimization')
parser.add_argument("--out", type=int, choices=[0,1], required=False, default=1, dest='out', help='Flag for outputting: 1 for observation and fit output')
parser.add_argument("--maxes", nargs='+', type=int, required=False, default=[2], dest='maxes', help='Indices of concentration maxima to include in error')
parser.add_argument("--yields", nargs='+', type=int, required=False, default=[], dest='yields', help='Indices of concentration yields to include in error')

args = parser.parse_args()

start=timeit.default_timer()
gas = ct.Solution(args.mechanism)
reactions=gas.reactions()
Npoints = args.Npoints
filebase = args.filebase
ns=gas.n_total_species
nr=gas.n_reactions
expfile = open(args.experiments, 'r')
measure_ind = args.measure
num_remove = args.remove
num_exclude = args.retain
maxes=args.maxes
yields=args.yields
seed = 1000*1000*num_exclude + 1000*num_remove + args.seed
ktol = args.ktol
outflag = args.out
filebase = args.filebase

if maxes == [-1]:
    maxes=[]
if yields == [-1]:
    yields=[]

tmaxes=[]
temperatures=[]
pressures=[]
initials=[]
for line in expfile:
    vals=line.split()
    tmaxes.append(float(vals[0]))
    temperatures.append(float(vals[1]))
    pressures.append(float(vals[2])*ct.one_atm)
    multiindex=np.zeros(ns)
    ics=[float(vals[m]) for m in range(3,len(vals))]
    for n in range(0,len(ics),2):
        multiindex[int(ics[n])]=ics[n+1]
    initials.append(multiindex)
nmeasure=len(temperatures)

observations=[]

if (not (os.path.exists('%sms.dat'%filebase) and os.path.exists('%sms.dat'%filebase))) and (outflag==1):
    print("running with sensitivities")
    for i in range(0,nmeasure):
        print("\n",i)
        observations.append(runsim(tmaxes[i], temperatures[i], pressures[i], initials[i]))
    sensitivities=[]
    for i in range(0,nmeasure):
        sensitivities.append(observations[i].sens)
    mi, ms=maxsens(sensitivities)

    np.savetxt('%smi.dat'%filebase,mi,fmt='%i')
    np.savetxt('%sms.dat'%filebase,ms)

    for n in range(0,nmeasure):
        np.save('%sobssens_%i.npy'%(filebase,n),np.array(observations[n].sens))
else:
    if (not (os.path.exists('%smi.dat'%filebase) and os.path.exists('%sms.dat'%filebase))):
        print("Data files missing; run with sensitivity!")
        quit()
    for i in range(0,nmeasure):
        observations.append(runsim_nosense(tmaxes[i], temperatures[i], pressures[i], initials[i]))
    mi=np.loadtxt('%smi.dat'%filebase,dtype=int)
    ms=np.loadtxt('%sms.dat'%filebase)

if (outflag == 1):
    print(nmeasure, " experimentes")
    print("tmaxes: ", *tmaxes)
    print("temperatures: ", *temperatures)
    print("pressures: ", *pressures)
    print("Maximum sensitivity norm for exact model range: %f - %f"%(np.min(ms),np.max(ms)))
    for i in mi[:num_exclude]:
        print(i, reactions[i].equation, ms[i])
    sys.stdout.flush()

    f=open('%sspecies.dat'%(filebase),'w')
    print(*(gas.species_names), sep=' ', file=f)
    f.close()
    f=open('%sreactants.dat'%(filebase),'w')
    for reaction in gas.reactions():
        for spec in reaction.reactants.keys():
            print(spec,reaction.reactants[spec], file=f, end=' ')
        print('', file=f)
    f.close()
    f=open('%sproducts.dat'%(filebase),'w')
    for reaction in gas.reactions():
        for spec in reaction.products.keys():
            print(spec,reaction.products[spec], file=f, end=' ')
        print('', file=f)
    f.close()
    for n in range(0,nmeasure):
        np.save("%stimes_%i.npy"%(filebase,n),observations[n].t)
        np.save("%stemperatures_%i.npy"%(filebase,n),observations[n].T)
        np.save("%spressures_%i.npy"%(filebase,n),observations[n].P)
        np.save("%smoles_%i.npy"%(filebase,n),observations[n].X)
        np.save("%srates_%i.npy"%(filebase,n),observations[n].net_rates_of_progress)

sensnorm = 0
removed=remove( seed, num_remove, mi[:num_exclude] )
for ind in removed:
    if(outflag == 1):
        print( ind, reactions[ind].equation, ms[ind])
    sensnorm += ms[ind]

if(gas.reaction_type(measure_ind) == 4):
    k0=gas.reactions()[measure_ind].low_rate.pre_exponential_factor
else:
    k0=gas.reactions()[measure_ind].rate.pre_exponential_factor
try:
    xa, xb, xc, fa, fb, fc, nf=op.bracket(residual, xa=np.log10(0.5), xb=np.log10(2.0), args=(k0, observations, measure_ind, tmaxes, temperatures, pressures, initials, maxes, yields), grow_limit=1.5)
    brack=(xa, xb, xc)
    if(outflag==1):
        print("bracket found in %d calls: (%f %f %f)"%(nf, xa, xb, xc))
    result=op.minimize_scalar(residual, args=(k0, observations, measure_ind, tmaxes, temperatures, pressures, initials, maxes, yields), method='brent', bracket=brack, options={'xtol': ktol})
except Exception as error:
    print('failed')
    print(error)
    f=open('%snorms.dat'%filebase,'w')
    print(measure_ind, seed-(1000*1000*num_exclude + 1000*num_remove), num_remove, num_exclude, sensnorm, "failed", "failed", *removed, sep=' ', file=f)
    f.close()
    stop=timeit.default_timer()
    print ('runtime: %f'%(stop-start))
    sys.stdout.flush()

    exit()


if(result.success):
    k=k0*10**(result.x)
    f=open('%snorms.dat'%filebase,'w')
    print(measure_ind, seed-(1000*1000*num_exclude + 1000*num_remove), num_remove, num_exclude, sensnorm, k/k0, result.fun, *removed, sep=' ', file=f)

    f.close()

    if (outflag == 1):
        print(result)
        sys.stdout.flush()
        reactions = gas.reaction_equations()
        if(gas.reaction_type(measure_ind) == 1):
            newrate=ct.ElementaryReaction(gas.reactions()[measure_ind].reactants,gas.reactions()[measure_ind].products)
            newrate.rate=ct.Arrhenius(k, gas.reactions()[measure_ind].rate.temperature_exponent, gas.reactions()[measure_ind].rate.activation_energy)
            gas.modify_reaction(measure_ind, newrate)
        if (gas.reaction_type(measure_ind) == 2):
            newrate=ct.ThreeBodyReaction(reactants=gas.reactions()[measure_ind].reactants,products=gas.reactions()[measure_ind].products)
            newrate.efficiencies = gas.reactions()[measure_ind].efficiencies
            newrate.rate=ct.Arrhenius(k, gas.reactions()[measure_ind].rate.temperature_exponent, gas.reactions()[measure_ind].rate.activation_energy)
            gas.modify_reaction(measure_ind, newrate)
        if (gas.reaction_type(measure_ind) == 4):
            newrate=ct.FalloffReaction(reactants=gas.reactions()[measure_ind].reactants,products=gas.reactions()[measure_ind].products)
            newrate.efficiencies = gas.reactions()[measure_ind].efficiencies
            newrate.falloff = gas.reactions()[measure_ind].falloff
            newrate.low_rate=ct.Arrhenius(k, gas.reactions()[measure_ind].low_rate.temperature_exponent, gas.reactions()[measure_ind].low_rate.temperature_exponent, gas.reactions()[measure_ind].low_rate.activation_energy)
            newrate.high_rate=gas.reactions()[measure_ind].high_rate
            gas.modify_reaction(measure_ind, newrate)

        print(measure_ind, reactions[measure_ind])
        for n in range(0,nmeasure):
            states = runsim_nosense ( tmaxes[n], temperatures[n], pressures[n], initials[n] )
            np.save("%sfit_%i.npy"%(filebase,n),states.X)

else:
    f=open('%snorms.dat'%filebase,'w')
    print(measure_ind, seed-(1000*1000*num_exclude + 1000*num_remove), num_remove, num_exclude, sensnorm, "failed", "failed", *removed, sep=' ', file=f)

    f.close()

stop=timeit.default_timer()
print ('runtime: %f'%(stop-start))
sys.stdout.flush()
