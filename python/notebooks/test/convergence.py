# -*- coding: utf-8 -*-
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
 
from plotting import *
import operator






def plot_transmittance_error():
	outpath = "c:/projects/visus/data/01_transmittance_study"
	dataset = "manix"
	methods = ["deltatracking", "ratiotracking", "residualratiotrackingmin", "residualratiotrackingavg", "residualratiotrackingmax", "stratifiedsampling", "randomoffset" ]


	numSamples_list = np.linspace(0, 10000, 100)
	plotSamples((numSamples_list, [1.0/sqrt(n) for n in numSamples_list]), "1/sqrt(N)", plot='loglog', linestyle = "--", color = "lightgray" )

	for method in methods:
		(numSamples_list, error_list) = loadSamplesFromFile( "{}/{}_{}_error.txt".format(outpath, dataset, method) )
		plotSamples((numSamples_list, error_list), "{}".format(method), plot='loglog')


	plt.xlabel(r'number of samples', fontsize=12)
	plt.ylabel(r'RMS', fontsize=12)
	plt.title("convergence")

	plt.legend(loc="best")
	plt.show()


def plot_transmittance_error_over_time():
	outpath = "c:/projects/visus/data/01_transmittance_study"
	dataset = "manix"
	methods = ["deltatracking", "ratiotracking", "residualratiotrackingmin", "residualratiotrackingavg", "residualratiotrackingmax", "stratifiedsampling", "randomoffset" ]


	for method in methods:
		(time_list, error_list) = loadSamplesFromFile( "{}/{}_{}_error_over_time.txt".format(outpath, dataset, method) )
		plotSamples((time_list, error_list), "{}".format(method), plot='logy')


	plt.xlabel(r'time(s)', fontsize=12)
	plt.ylabel(r'RMS', fontsize=12)
	plt.title("convergence")

	plt.legend(loc="best")


	plt.show()



	'''
	(numSamples_list, error_raymarching) = loadSamplesFromFile( "c:/projects/msc/data/lineintegration_error_raymarching" )
	(numSamples_list, error_deltatracking) = loadSamplesFromFile( "c:/projects/msc/data/lineintegration_error_deltatracking" )
	#(numSamples_list, error_splitteddeltatracking) = loadSamplesFromFile( "c:/projects/msc/data/lineintegration_error_splitteddeltatracking" )
	(numSamples_list, error_ratiotracking) = loadSamplesFromFile( "c:/projects/msc/data/lineintegration_error_ratiotracking" )
	(numSamples_list, error_residualratiotracking) = loadSamplesFromFile( "c:/projects/msc/data/lineintegration_error_residualratiotracking" )
	(numSamples_list, error_stratifiedsampling) = loadSamplesFromFile( "c:/projects/msc/data/lineintegration_error_stratifiedsampling" )
	

	plotSamples((numSamples_list, error_raymarching), "raymarching", plot='logy')
	plotSamples((numSamples_list, error_deltatracking), "delta tracking", plot='logy')
	#plotSamples((numSamples_list, error_splitteddeltatracking), "splitted delta tracking", plot='logy')
	plotSamples((numSamples_list, error_ratiotracking), "ratio tracking", plot='logy')
	plotSamples((numSamples_list, error_residualratiotracking), "residual ratio tracking", plot='logy')
	plotSamples((numSamples_list, error_stratifiedsampling), "stratified sampling", plot='logy')
	plotSamples((numSamples_list, error_randomoffset), "random offset", plot='logy')
	#plotSamples( (x_list, y_list), "test", color='r' )

	#(t_list, pdf_list) = loadSamplesFromFile( "/disney/users/dkoerner/temp/svpt_renders/exit_sampling/subdiv/samples_pdf" );
	#plotSamples( (t_list, pdf_list), "test", color='r', marker='.', linestyle=' ' )

	plt.xlabel(r'number of samples', fontsize=12)
	plt.ylabel(r'RMS', fontsize=12)
	plt.title("convergence")

	plt.legend(loc="best")
	plt.show()


	(numSamples_list, time_raymarching) = loadSamplesFromFile( "c:/projects/msc/data/lineintegration_time_raymarching" )
	(numSamples_list, time_deltatracking) = loadSamplesFromFile( "c:/projects/msc/data/lineintegration_time_deltatracking" )
	#(numSamples_list, time_splitteddeltatracking) = loadSamplesFromFile( "c:/projects/msc/data/lineintegrationtime_splitteddeltatracking" )
	(numSamples_list, time_ratiotracking) = loadSamplesFromFile( "c:/projects/msc/data/lineintegration_time_ratiotracking" )
	(numSamples_list, time_residualratiotracking) = loadSamplesFromFile( "c:/projects/msc/data/lineintegration_time_residualratiotracking" )
	(numSamples_list, time_stratifiedsampling) = loadSamplesFromFile( "c:/projects/msc/data/lineintegration_time_stratifiedsampling" )
	(numSamples_list, time_randomoffset) = loadSamplesFromFile( "c:/projects/msc/data/lineintegration_time_randomoffset" )

	plotSamples((numSamples_list, time_raymarching), "raymarching", plot='logy')
	plotSamples((numSamples_list, time_deltatracking), "delta tracking", plot='logy')
	#plotSamples((numSamples_list, time_splitteddeltatracking), "splitted delta tracking", plot='logy')
	plotSamples((numSamples_list, time_ratiotracking), "ratio tracking", plot='logy')
	plotSamples((numSamples_list, time_residualratiotracking), "residual ratio tracking", plot='logy')
	plotSamples((numSamples_list, time_stratifiedsampling), "stratified sampling", plot='logy')
	plotSamples((numSamples_list, time_randomoffset), "random offset", plot='logy')
	#plotSamples( (x_list, y_list), "test", color='r' )

	plt.xlabel(r'number of samples', fontsize=12)
	plt.ylabel(r'time (seconds)', fontsize=12)
	plt.title("performance")
	'''



def plot_transmittance_variance_over_time():
	outpath = "c:/projects/msc/data/bias_test2"
	dataset = "sin"
	methods = ["ro", "rt" ]


	for method in methods:
		(time_list, error_list) = loadSamplesFromFile( "{}/lineintegration_{}_{}_variance_over_time.txt".format(outpath, dataset, method) )
		plotSamples((time_list, error_list), "{}".format(method), plot='logy')
	(time_list, error_list) = loadSamplesFromFile( "{}/lineintegration_{}_{}_variance_over_time_cv.txt".format(outpath, dataset, "ro") )
	plotSamples((time_list, error_list), "{}_cv".format("ro"), plot='logy')


	plt.xlabel(r'time(ms)', fontsize=12)
	plt.ylabel(r'variance', fontsize=12)
	#plt.title("convergence")

	plt.legend(loc="best")

	pp = PdfPages('sin_variance_over_time.pdf')
	plt.savefig(pp, format='pdf', facecolor='white', bbox_inches='tight')
	pp.close() 

	plt.show()

def plot_transmittance_error_over_time_allrays():
	outpath = "c:/projects/visus/data/01_transmittance_study"
	#datasets = [ "foraminifera", "chameleon", "nebulae", "manix", "artifix"]
	datasets = [ "manix"]
	for dataset in datasets:
		#dataset = "nebulae"
		methods = ["deltatracking", "ratiotracking", "residualratiotrackingmin", "residualratiotrackingavg", "residualratiotrackingmax", "stratifiedsampling", "randomoffset", "randomoffsettestns1", "randomoffsettestns2", "randomoffsettestns3" ]


		for method in methods:
			(time_list, error_list) = loadSamplesFromFile( "{}/{}_{}_error_over_time_allrays.txt".format(outpath, dataset, method) )
			if method == "randomoffsettest":
				plotSamples((time_list, error_list), "{}".format(method), plot='logy', marker='.')
			else:
				plotSamples((time_list, error_list), "{}".format(method), plot='logy')


		plt.xlabel(r'time(s)', fontsize=12)
		plt.ylabel(r'RMS', fontsize=12)
		plt.title("convergence | {}".format(dataset))

		plt.legend(loc="best")

		pp = PdfPages('manix_error_over_time_all_rays.pdf')
		plt.savefig(pp, format='pdf', facecolor='white', bbox_inches='tight')
		pp.close() 

		plt.show()

def plot_transmittance_error_over_samples_allrays():
	outpath = "c:/projects/visus/data/01_transmittance_study"
	#datasets = [ "foraminifera", "chameleon", "nebulae", "manix", "artifix"]
	datasets = [ "manix"]
	for dataset in datasets:
		#dataset = "nebulae"
		#methods = ["deltatracking", "ratiotracking", "residualratiotrackingmin", "residualratiotrackingavg", "residualratiotrackingmax", "stratifiedsampling", "randomoffset" ]
		methods = ["deltatracking", "ratiotracking", "randomoffset", "randomoffsettest" ]


		for method in methods:
			(samples_list, error_list) = loadSamplesFromFile( "{}/{}_{}_error_over_samples_allrays.txt".format(outpath, dataset, method) )
			plotSamples((samples_list, error_list), "{}".format(method), plot='logy')


		plt.xlabel(r'#samples', fontsize=12)
		plt.ylabel(r'RMS', fontsize=12)
		plt.title("convergence | {}".format(dataset))

		plt.legend(loc="best")

		#pp = PdfPages('manix_error_over_time_all_rays.pdf')
		#plt.savefig(pp, format='pdf', facecolor='white', bbox_inches='tight')
		#pp.close() 

		plt.show()

def plot_ray():

	#name = "transmittance"

	name = "od"
	(index_list, sample_list) = loadSamplesFromFile("c:/projects/visus/data/01_transmittance_study/manix_analysis_ray_{}".format(name))
	plotSamples((index_list, sample_list), name)

	#for i in range(10):
	numSegments = 10
	dt = 7.95153
	for i in range(numSegments+1):
		plt.axvline(i*dt, color="lightgray")

	plt.legend()
	plt.show()


def computeT( extinction_list, extinction_max, muc = 0.0 ):
	T = 1.0
	for extinction in extinction_list:
		T = T*(1-(extinction-muc)/extinction_max)
	return T



def plot_ratiotracking_test():

	(step_index, extinction_list) = loadSamplesFromFile("c:/projects/visus/data/01_transmittance_study_test/test")
	extinction_max = 0.81868
	numSteps = len(step_index)
	T = 1.0
	T_list = []
	for i in range(numSteps):
		T = computeT(extinction_list[:i+1], extinction_max)
		T_list.append( T )

	muc_list = np.linspace(0, extinction_max, 30)
	test_list = []
	for muc in muc_list:
		T = computeT(extinction_list[:i+1], extinction_max, muc*extinction_max)
		test_list.append(T)

	#plotSamples((step_index, T_list), "")
	plotSamples((muc_list, test_list), "")


	#plotSamples((muc_list, [exp(x*7) for x in muc_list]), "")

	plt.show()



def gauss( x, mu, stddev ):
	t = (x-mu)/stddev
	return 1.0/(stddev*sqrt(2*pi))*exp(-0.5*t*t)


def plot_1( basepath, title ):

	(x_list_extinction, extinction_groundtruth_list) = loadSamplesFromFile("{0}/debug_extinction_1.txt".format(basepath))
	(x_list_svd, extinction_svd_list) = loadSamplesFromFile("{0}/debug_extinction_3.txt".format(basepath))
	(temp, extinction_svd_segments_list) = loadSamplesFromFile("{0}/debug_extinction_3_segments.txt".format(basepath))
	(x_list_segmentidx, segmentidx_list) = loadSamplesFromFile("{0}/debug_segmentidx_3.txt".format(basepath))


	(x_list_extinction_segment, extinction_groundtruth_list_segment) = loadSamplesFromFile("{0}/debug_segment_extinction_1.txt".format(basepath))
	(x_list_extinction_segment, extinction_svd_list_segment) = loadSamplesFromFile("{0}/debug_segment_extinction_3.txt".format(basepath))
	plotSamples( (x_list_extinction_segment, extinction_groundtruth_list_segment), "groundtruth_segment")
	plotSamples( (x_list_extinction_segment, extinction_svd_list_segment), "svd_segment")


	#ax = plt.gca()
	#ax.grid(True, which='both')

	#for s in extinction_svd_segments_list:
	#	ax.axvline(x=s, color='k')

	#ax.axhline(y=0, color='k')
	#ax.axvline(x=0, color='k')
	plt.title(title)
	plt.legend()
	plt.show()

if __name__ == "__main__":

	#plot_transmittance_error_over_samples()
	#plot_transmittance_error_over_time()
	#plot_transmittance_error_over_time_allrays()
	#plot_transmittance_error_over_samples_allrays()
	#plot_transmittance_variance_over_time()

	# fitting ----
	#(x_list, y_list) = loadSamplesFromFile("c:/projects/visus/data/noisereduction/fitting_gradient/C60_input_distr.txt")
	#y_list = y_list/np.max(y_list)
	#plt.plot(x_list, y_list, "o")



	#x_list = np.linspace(0, 2, 200)
	#mu = 0.0359
	#mu = 0.0
	#sigma = 0.0744
	#gauss_max = gauss(mu, mu, sigma)
	#y_list = [gauss(x, mu, sigma)/gauss_max for x in x_list]
	#plotSamples( (x_list, y_list), "test", color='r' )

	#plt.show()



	#plot_ray()

	#plot_ratiotracking_test()


	#plot_1("c:/projects/visus/doc/meetings/10_fit_svd_bias/weight_from_optimization", "weight=optimization") # weighting optimization (best match of "mass")
	#plot_1("c:/projects/visus/doc/meetings/10_fit_svd_bias/centroid_distance_threshold_01", "weight=max_extinction | samples with distance < 0.01 | segment length = 0.024")
	#plot_1("c:/projects/visus/doc/meetings/10_fit_svd_bias/centroid_distance_threshold_025", "weight=max_extinction | samples with distance < 0.025 | segment length = 0.024")

	#plot_1("c:/projects/visus/data/noisereduction/fitting3_svd", "")

	basepath = "c:/projects/epfl/epfl17/cpp/build-adrrs-Desktop_Qt_5_4_1_MSVC2013_OpenGL_64bit-Release";
	order = 6
	for l in range(order):
		for m in range(0, l+1):
			print("l={} m={}".format(l,m))
			if not( l == 5):
				continue;
			#if not(m == 0 and l==2):
			#	continue;
			(x_list, y_list) = loadSamplesFromFile("{}/test_P_{}_{}.txt".format(basepath, l, m))
			(x_list, y2_list) = loadSamplesFromFile("{}/test_P2_{}_{}.txt".format(basepath, l, m))
			ax, = plotSamples( (x_list, y_list), label="{} {}".format(l,m))
			plotSamples( (x_list, y2_list), linestyle=" ", marker=".", color = ax.get_color())

	#ax = plt.gca()
	#ax.grid(True, which='both')

	#for s in extinction_svd_segments_list:
	#	ax.axvline(x=s, color='k')

	#ax.axhline(y=0, color='k')
	#ax.axvline(x=0, color='k')
	#plt.title(title)
	plt.legend(loc='best')
	plt.show()

	sys.exit(0)