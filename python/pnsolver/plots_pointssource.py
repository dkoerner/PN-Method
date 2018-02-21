import os,sys,inspect
import random

import renderer
import util

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

sigma_t = 8.0
albedo = 0.9
size = 2.0
pointsource_center = np.array([1.0, 1.0, 1.0])

def grosjean(r, sigma_t, albedo):
    sigma_a = (1.0-albedo)*sigma_t
    sigma_s = albedo*sigma_t
    ss = np.exp(-sigma_t*r)/(4.0*np.pi*r*r)
    tmp = sigma_a+sigma_s
    D = (2.0*sigma_a+sigma_s)/(3.0*tmp*tmp)
    ms = 3.0*sigma_s*sigma_t*np.exp(-np.sqrt(sigma_a/D)*r)/(4.0*np.pi*(2.0*sigma_a+sigma_s)*r)
    return ss+ms
    #return ms

def grosjean_ss(r, sigma_t, albedo):
    sigma_a = (1.0-albedo)*sigma_t
    sigma_s = albedo*sigma_t
    ss = np.exp(-sigma_t*r)/(4.0*np.pi*r*r)
    tmp = sigma_a+sigma_s
    D = (2.0*sigma_a+sigma_s)/(3.0*tmp*tmp)
    ms = 3.0*sigma_s*sigma_t*np.exp(-np.sqrt(sigma_a/D)*r)/(4.0*np.pi*(2.0*sigma_a+sigma_s)*r)
    #return ss+ms
    return ss

def cda_org( r, sigma_t, albedo ):
    sigma_a = (1.0-albedo)*sigma_t
    sigma_tr = np.sqrt(3.0*sigma_a*sigma_t)
    return 3.0*sigma_t*np.exp(-sigma_tr*r)/(4.0*np.pi*r)


pns_cda = renderer.load_pnsolution( "C:/projects/epfl/epfl17/python/pnsolver/results/pointsource/pointsource_cda5.pns" )
pns_cda_old = renderer.load_pnsolution( "C:/projects/epfl/epfl17/python/pnsolver/results/pointsource/pointsource_cda_old.pns" )
pns_cda_new = renderer.load_pnsolution( "C:/projects/epfl/epfl17/python/pnsolver/results/pointsource/pointsource_cda_new2.pns" )
pns_fld_new = renderer.load_pnsolution( "C:/projects/epfl/epfl17/python/pnsolver/results/pointsource/pointsource_fld_new.pns" )

pns_p1_old = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/pointsource/pointsource_p1_old.pns" )
pns_p1_new = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/pointsource/pointsource_p1_new.pns" )

pns_p2 = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/pointsource/pointsource_p2_new.pns" )
pns_p3 = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/pointsource/pointsource_p3_new.pns" )
pns_p4 = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/pointsource/pointsource_p4_new.pns" )
pns_p5 = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/pointsource/pointsource_p5_new.pns" )

pns_p1_ms = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/pointsource/pointsource_p1_new_ms.pns" )

#pns_p5 = renderer.load_pnsolution( "c:/projects/epfl/epfl17/python/pnsolver/results/pointsource/pointsource_p5_4_ms.pns" )

r_list = np.linspace(1.0e-2, size*0.5, 100)
fluence_grosjean = [ grosjean(r, sigma_t, albedo) for r in r_list ]
fluence_cda = [ cda_org(r, sigma_t, albedo) for r in r_list ]

fluence_cda_old = np.array([pns_cda_old.evalCoefficient(pointsource_center+np.array([r, 0.0, 0.0]), 0) for r in r_list])
fluence_cda_new = np.array([pns_cda_new.evalCoefficient(pointsource_center+np.array([r, 0.0, 0.0]), 0) for r in r_list])
fluence_fld_new = np.array([pns_fld_new.evalCoefficient(pointsource_center+np.array([r, 0.0, 0.0]), 0) for r in r_list])

#fluence_p1_numerical = np.array([pns_p1.evalCoefficient(pointsource_center+np.array([r, 0.0, 0.0]), 0) for r in r_list])
fluence_p1_old = np.array([pns_p1_old.evalCoefficient(pointsource_center+np.array([r, 0.0, 0.0]), 0) for r in r_list])
fluence_p1_new = np.array([pns_p1_new.evalCoefficient(pointsource_center+np.array([r, 0.0, 0.0]), 0) for r in r_list])
fluence_p2 = np.array([pns_p2.evalCoefficient(pointsource_center+np.array([r, 0.0, 0.0]), 0) for r in r_list])
fluence_p3 = np.array([pns_p3.evalCoefficient(pointsource_center+np.array([r, 0.0, 0.0]), 0) for r in r_list])
fluence_p4 = np.array([pns_p4.evalCoefficient(pointsource_center+np.array([r, 0.0, 0.0]), 0) for r in r_list])
fluence_p5 = np.array([pns_p5.evalCoefficient(pointsource_center+np.array([r, 0.0, 0.0]), 0) for r in r_list])


t = np.sqrt(4.0*np.pi)
fluence_p1_ms = np.array([grosjean_ss(r, sigma_t, albedo)+pns_p1_ms.evalCoefficient(pointsource_center+np.array([r, 0.0, 0.0]), 0) for r in r_list])
#fluence_p5_numerical = np.array([pns_p5.evalCoefficient(pointsource_center+np.array([r, 0.0, 0.0]), 0) for r in r_list])





'''
fig = plt.figure(figsize=(8,8));
ax = fig.add_subplot(111)
plt.loglog(r_list, fluence_grosjean, label="Grosjean", color="k")
#plt.loglog(r_list, fluence_cda, label="CDA", color="r")
#plt.loglog(r_list, fluence_cda_old, label="CDA old", color="r")
#plt.loglog(r_list, fluence_cda_new, label="CDA new", color="r")
#plt.loglog(r_list, fluence_cda_numerical, label="CDA numerical", marker=".", linestyle=" " )

#plt.loglog(r_list, fluence_p1_numerical, label="P1", marker=" ", linestyle="-", color=".95" )
plt.loglog(r_list, fluence_p1_new, label="P1", marker=" ", linestyle="-", color=(0.0,.95,0.0) )
#plt.loglog(r_list, fluence_p1_numerical_check, label="P1 check", marker=".", linestyle=" ", color=".95" )
plt.loglog(r_list, fluence_p2, label="P2", marker=" ", linestyle="-", color=(0.0,.8,0.0) )
plt.loglog(r_list, fluence_p3, label="P3", marker=" ", linestyle="-", color=(0.0,.6,0.0) )
plt.loglog(r_list, fluence_p4, label="P4", marker=" ", linestyle="-", color=(0.0,.4,0.0) )
plt.loglog(r_list, fluence_p5, label="P5", marker=" ", linestyle="-", color=(0.0,.2,0.0) )
plt.xlabel('radius')
plt.ylabel('fluence')
plt.legend()
plt.savefig('pointsource_pn.pdf', bbox_inches='tight')
plt.show()
'''






fig = plt.figure(figsize=(8,8));
ax = fig.add_subplot(111)
plt.loglog(r_list, fluence_grosjean, label="Grosjean", color="k")
plt.loglog(r_list, fluence_cda_new, label="CDA", color="r")
plt.loglog(r_list, fluence_fld_new, label="FLD", color="b")
#plt.loglog(r_list, fluence_cda_old, label="CDA old", color="r")
#plt.loglog(r_list, fluence_cda_new, label="CDA new", color="r")
#plt.loglog(r_list, fluence_cda_numerical, label="CDA numerical", marker=".", linestyle=" " )

#plt.loglog(r_list, fluence_p1_numerical, label="P1", marker=" ", linestyle="-", color=".95" )
plt.loglog(r_list, fluence_p1_new, label="P1", marker=".", linestyle=" ", color="r" )
#plt.loglog(r_list, fluence_p1_numerical_check, label="P1 check", marker=".", linestyle=" ", color=".95" )
#plt.loglog(r_list, fluence_p2, label="P2", marker=" ", linestyle="-", color=".8" )
#plt.loglog(r_list, fluence_p3, label="P3", marker=" ", linestyle="-", color=".6" )
#plt.loglog(r_list, fluence_p4, label="P4", marker=" ", linestyle="-", color=".4" )
plt.loglog(r_list, fluence_p5, label="P5", marker=" ", linestyle="-", color="g" )
#plt.loglog(r_list, fluence_p1_ms, label="P1 ms", marker=".", linestyle=" ", color="g" )
plt.xlabel('radius')
plt.ylabel('fluence')
plt.legend()
#plt.savefig('pointsource_p5.pdf', bbox_inches='tight')
plt.show()