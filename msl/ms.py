from random import sample
from random import shuffle

from math import floor, ceil, sqrt

import networkx as nx
from networkx.algorithms.approximation import clique

from blincodes import matrix, vector
from blincodes.codes import rm, tools

class MinderShokrollahi():
	''' Class which contains functions to perform one step of Minder-Shokrollahi attack '''

	def __init__(self, r, m):
		self.r = r
		self.m = m
		self.d = 2**(self.m - self.r)

		self.pair_ctrs = []
		for i in range(2**self.m):
			self.pair_ctrs.append([0 for j in range(2**self.m)])	

		self.min_wt_cws = []

	def min_weight_sample(self, generator):
		''' Finds minimum weight codeword '''	

		for row in generator:
			if row.support not in self.min_wt_cws and row.hamming_weight == self.d:
				self.min_wt_cws.append(row.support)
				return row.support

		while True:
			cols_sample = sample(range(generator.ncolumns), generator.nrows)
			gaussed_g = generator.gaussian_elimination(cols_sample)

			for row in gaussed_g:
				if row.support not in self.min_wt_cws and row.hamming_weight == self.d:
					self.min_wt_cws.append(row.support)
					return row.support


	def get_cliques(self, G):

		clique_size = self.d

		result = []
		not_found = True

		while (True):

			cliques = nx.find_cliques(G)

			for clique in cliques:
				
				if (not_found):
					this_clique_size = len(clique)

					if this_clique_size >= clique_size:
						if not this_clique_size % clique_size:

							G.remove_nodes_from(clique)
							not_found = False

							for i in range(int(this_clique_size / clique_size)):

								result.append(clique[i*clique_size: (i + 1)*clique_size])

			if len(result) == 2**self.r - 1:
				return result

			if not_found:
				return []
			else:
				not_found = True

	def decompose_inner_sets(self, generator):

		cw_num = 80 * 2**(self.m - 5)
		threshold = 100

		eps = (1 - 2**(-self.r)) * sqrt(1 - 2**(self.r - 1) / self.d)
		min_weight = self.d
		max_weight = floor(2 * self.d * eps)

		low_wt_cws = []

		while len(low_wt_cws) < cw_num:
			cols_sample = sample(range(generator.ncolumns), generator.nrows)
			new_generator = matrix.nonsingular(generator.nrows) * generator
			gaussed_g = new_generator.gaussian_elimination(cols_sample)

			for row in gaussed_g:
				wt = row.hamming_weight
				if row.support not in low_wt_cws and wt >= min_weight and wt <= max_weight:
					low_wt_cws.append(row.support)

		for i in range(2**self.m):
			for j in range(2**self.m):
				self.pair_ctrs[i][j] = 0


		for codeword in low_wt_cws:
			for i in range(2**self.m):
				for j in range(i + 1, 2**self.m):
					if i in codeword and j in codeword:
						self.pair_ctrs[i][j] += 1


		G = nx.Graph()
		G.add_nodes_from(range(2**self.m))
		cliques = []

		iteration = 0
		max_iterations = 1000

		while not len(cliques) and iteration < max_iterations:

			for i in range(2**self.m):
				for j in range(i + 1, 2**self.m):
					if (self.pair_ctrs[i][j] * (iteration + 1)) >= threshold:
						G.add_edge(i, j)

			cliques = self.get_cliques(G)

			iteration += 1

		return cliques

	def binomial_coef(self, n, k):
		numerator_prod = 1
		denominator_prod = 1
		for i in range(k):
			numerator_prod *= n - i
			denominator_prod *= i + 1

		return numerator_prod // denominator_prod

	def attack(self, generator):
		rm_subcode_basis = matrix.Matrix()
		rm_subcode_dim = 0
		for i in range(self.r):
			rm_subcode_dim += self.binomial_coef(self.m, i)

		while rm_subcode_basis.nrows < rm_subcode_dim: 
			min_weight_codeword = self.min_weight_sample(generator)

			shortened_g = tools.truncate(generator, min_weight_codeword)

			inner_sets = self.decompose_inner_sets(shortened_g)

			f_vecs = [vector.from_support(2**self.m, min_weight_codeword + s) for s in inner_sets]

			rm_subcode_basis = tools.union(rm_subcode_basis, matrix.from_vectors(f_vecs))


		return rm_subcode_basis
