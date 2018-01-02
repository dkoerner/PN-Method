# -*- coding: utf-8 -*-


import math
from matplotlib.pylab import *
import matplotlib.pyplot  as pyplot
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
import numpy as np
from functools import partial



def interpolate( x_list, y_list, x ):
	if x>=x_list[-1]:
		return y_list[-1]
	if x<=x_list[0]:
		return y_list[0]

	numSamples = len(x_list)
	start = 0
	for i in range(numSamples-1):
		if x < x_list[i+1]:
			start = i
			break
	t = (x - x_list[start])/(x_list[start+1]-x_list[start])
	return t*y_list[start+1] + (1.0-t)*y_list[start]

def plot_sample_distribution( samples, bin = 0.01, normalizeTo = 1.0, label = "distribution", color=None ):
	x = samples
	# find min, max range
	maxt = -999999.0
	mint =  999999.0
	for i in range(len(x)):
		mint = min(mint, x[i])
		maxt = max(maxt, x[i])

	# num segments
	d = bin
	numSegments = (maxt-mint)/d
	
	segments = []
	dist = []
	# maximum number a single segment has in the distribution
	# used to renormalize distribution values to 0-1
	maxSamplesPerSegment = 0
	for s in samples:
		segment = int((s-mint)/d)

		while segment >= len(dist):
			dist.append(0)
		dist[segment] += 1
		
		if dist[segment] > maxSamplesPerSegment:
			maxSamplesPerSegment = dist[segment]

	# normalize
	for i in range(len(dist)):
		segments.append(mint + i*d)
		dist[i] = float(dist[i])/float(maxSamplesPerSegment)
		dist[i] *= normalizeTo

	#plt.plot(segments, dist, label=label, linewidth=2)
	plt.bar(segments, dist, bin, color=color)

def plotSamples(samples, label="", plot=None, redistributeSamples=False, linewidth=1.0, color = None, scale = 1.0, scale_x = 1.0, marker= 'None', linestyle = '-', lowerbound = None, markersize=10 ):

	x = samples[0]
	y = samples[1]
	
	#scale y
	for i in range( len(y) ):
		x[i] = x[i]*scale_x
		if y[i] != None:
			y[i] = y[i]*scale


	if (plot == 'logx' or plot == 'loglog') and redistributeSamples==True:
		# redistribute x-axis samples in logspace such that samples are evenly distributed
		numSamples = 50
		x_log = np.logspace(-2, -0.29, num = numSamples, endpoint=True)
		y_log = []
		for i in range( len(x_log) ):
			y_log.append( interpolate( x, y, x_log[i] ) )
		
		x = x_log
		y = y_log
	
	if color != None:
		if plot == None:
			return plt.plot(x,y, label=label, linewidth=linewidth, color=color, marker=marker, linestyle=linestyle, markersize=markersize,markeredgecolor=color)
		elif plot == 'logy':
			return plt.semilogy(x,y, label=label, linewidth=linewidth, color=color, marker=marker, linestyle=linestyle, markersize=markersize,markeredgecolor=color)
		elif plot == 'logx':
			return plt.semilogx(x,y, label=label, linewidth=linewidth, color=color, marker=marker, linestyle=linestyle, markersize=markersize,markeredgecolor=color)
		elif plot == 'loglog':
			return plt.loglog(x,y, label=label, linewidth=linewidth, color=color, marker=marker, linestyle=linestyle, markersize=markersize,markeredgecolor=color) #markeredgewidth=0.0
	else:
		if plot == None:
			return plt.plot(x,y, label=label, linewidth=linewidth, marker=marker, linestyle=linestyle, markersize=markersize,markeredgecolor=color)
		elif plot == 'logy':
			return plt.semilogy(x,y, label=label, linewidth=linewidth, marker=marker, linestyle=linestyle, markersize=markersize,markeredgecolor=color)
		elif plot == 'logx':
			return plt.semilogx(x,y, label=label, linewidth=linewidth, marker=marker, linestyle=linestyle, markersize=markersize,markeredgecolor=color)
		elif plot == 'loglog':
			return plt.loglog(x,y, label=label, linewidth=linewidth, marker=marker, linestyle=linestyle, markersize=markersize,markeredgecolor=color)

def loadSamplesFromFile( filename, lowerbound = None ):
	print( "loading samples from {}".format(filename) )
	data = []
	f = open(filename,'r')
	fdata = f.read()
	f.close()

	rows = fdata.split('\n')
	x = []
	y = []
	for column in rows:
		columns = column.split(' ')
		if len(columns) == 1:
			if len(columns[0])>0:
				x.append(float(columns[0]))
		if len(columns) == 2:
			x.append(float(columns[0]))

			yyy = None
			try:
				yyy = float(columns[1])
				if lowerbound:
					if yyy < lowerbound:
						yyy = 0.0
			except Exception(e):
				pass
			
			y.append(yyy)
	return (np.array(x), np.array(y))
	
def loadSurfaceFromFile( filename ):
	data = []
	f = open(filename,'r')
	info = f.readline()
	uv = info.split(' ')
	if len(uv) == 2:
		u = int(uv[0])
		v = int(uv[1])
	else:
		return
	
	data = f.read()
	f.close()
	points = data.split('\n')

	# first line in the file gives the number of points in u and v on the 2d manifold of the surface
	x = []
	y = []
	z = []
	for point in points:
		xyz = point.split(' ')
		if len(xyz) == 3:
			x.append(float(xyz[0]))
			y.append(float(xyz[1]))
			z.append(float(xyz[2]))

	X = np.array(x).reshape(u,v)
	Y = np.array(y).reshape(u,v)
	Z = np.array(z).reshape(u,v)
	return (X,Y,Z)


	
def plotFile(filename, label="", plot=None, redistributeSamples=False, linewidth=1.0, color = None, scale = 1.0, scale_x = 1.0, marker= None, linestyle = '-', lowerbound = None, markersize=10 ):
	
	(x, y) = loadSamplesFromFile(filename, lowerbound)	
	plotSamples( (x, y), label, plot, redistributeSamples, linewidth, color, scale, scale_x, marker, linestyle, lowerbound, markersize )
	
	return (x, y)
	
	
def sample( fun, domain_a, domain_b, numSamples=10 ):
	x = []
	y = []
	for i in range(numSamples):
		x_value = (domain_a+(domain_b-domain_a)/float(numSamples-1)*float(i))
		x.append(x_value)
		y.append(fun(x_value))
	return (x, y)
	
def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)