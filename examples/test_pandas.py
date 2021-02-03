# @fileoverview Copyright (c) 2019-2020, Stefano Gualandi,
#               via Ferrata, 1, I-27100, Pavia, Italy
#
# @author stefano.gualandi@gmail.com (Stefano Gualandi)


# Library for reading input files
import pickle
import numpy

from KWD import *

# open a file, where you stored the pickled data
with open('Iteration_1_250.pickle', 'rb') as file:

	# read information from that file
	data = pickle.load(file)

	N = 2000

	# Scale coordinates
	X = ((data['i'] - min(data['i']))/125).to_numpy(dtype=int)[:N]
	Y = ((data['j'] - min(data['j']))/125).to_numpy(dtype=int)[:N]
	
	# Data histograms
	H1 = data['h1'].to_numpy(dtype=float)[:N]
	H2 = data['h2'].to_numpy(dtype=float)[:N]

	# Numebr of items
	n = len(X)

	# Define first histogram
	a = Histogram2D(n, X, Y, H1)

	# Define second histogram
	b = Histogram2D(n, X, Y, H2)

	# Define a solver and compute the distance 
	# Kantorovich-Wasserstein distance of order 1
	# with L_2 as ground distance
	s = Solver()

	print("approx: d(a,b) =", s.distance(a, b, 2), "- runtime:", s.runtime())
	print("approx: d(a,b) =", s.distance(a, b, 3), "- runtime:", s.runtime())
	print("approx: d(a,b) =", s.distance(a, b, 5), "- runtime:", s.runtime())

	print("exact:  d(a,b) =", s.dense(a, b), "- runtime:", s.runtime())