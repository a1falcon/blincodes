from random import shuffle
from math import ceil

from blincodes import matrix, vector
from blincodes.codes import rm

from msl import ms


class ChizhovBorodin:
	''' Class which contains functions to perform Chizhov-Borodin attack '''

	def __init__(self, r, m):
		self.r = r
		self.m = m

	def circle_dot_prod(self, generator1, generator2):

		product_list = [] 
		for row1 in generator1:
			for row2 in generator2:
				product_list.append(row1 & row2)


		prod_g = matrix.from_vectors(product_list)

		prod_g = prod_g.gaussian_elimination()

		res_g = []
		for row in prod_g:
			if row.hamming_weight:
				res_g.append(row)

		return matrix.from_vectors(res_g)


	def pubkey_gen(self):
		rm_generator = rm.generator(self.r, self.m)
		M = matrix.nonsingular(rm_generator.nrows)
		perm = [i for i in range(rm_generator.ncolumns)]
		shuffle(perm)
		P = matrix.permutation(perm)

		return M * rm_generator * P

	def extended_gcd(self, a, b):
		''' Finds x, y and d such that ax + by = d = gcd(a,b) '''
		r_prev, r_cur = a, b
		x_prev, x_cur = 1, 0
		y_prev, y_cur = 0, 1

		while r_cur:
			q = r_prev // r_cur
			r_prev, r_cur = r_cur, r_prev - q * r_cur
			x_prev, x_cur = x_cur, x_prev - q * x_cur
			y_prev, y_cur = y_cur, y_prev - q * y_cur

		return (x_prev, y_prev, r_prev)

	def positive_a_case(self, a, b, generator):
		q = ceil(-b / a)
		s = q * a + b

		if not s:

			rm_qr_orth = generator

			for i in range(0,q-1):
				rm_qr_orth = self.circle_dot_prod(rm_qr_orth, generator)

			rm_qr_orth = rm_qr_orth.orthogonal

			rm_qr_orth_mul = rm_qr_orth

			for i in range(0,a-1):
				rm_qr_orth_mul = self.circle_dot_prod(rm_qr_orth_mul, rm_qr_orth)

			return rm_qr_orth_mul

		elif q and s:

			rm_sr = generator

			for i in range(0,s-1):
				rm_sr = self.circle_dot_prod(rm_sr, generator)

			rm_qr_orth = generator

			for i in range(0,q-1):
				rm_qr_orth = self.circle_dot_prod(rm_qr_orth, generator)

			rm_qr_orth = rm_qr_orth.orthogonal

			rm_qr_orth_mul = rm_qr_orth

			for i in range(0,a-1):
				rm_qr_orth_mul = self.circle_dot_prod(rm_qr_orth_mul, rm_qr_orth)

			return self.circle_dot_prod(rm_qr_orth_mul, rm_sr)	

	def generate_rm_d(self, generator, a, b, d):

		if not a and b == 1:
			return generator

		if a > 0 and b < 0:
			return self.positive_a_case(a, b, generator)

		if a < 0 and b > 0:
			return self.positive_a_case(1-a, -b, generator).orthogonal

	def find_permutation(self, generator):

		onev = vector.from_support_supplement(2**m)

		a = generator.T.solve(onev)[1]

		removing_num = 0

		if len(a.support):
			removing_num = a.support[0]

		A_rows = [a]

		for i in range(0, m + 1):
			if i != removing_num:
				A_rows.append(a ^ vector.from_support(m + 1, [i]))

		ag = matrix.from_vectors(A_rows)*generator

		A_rows = ag[1:]
		
		return matrix.permutation([row.value for row in A_rows.T])


	def find_nonsingular(self, generator, g_mul_perm):

		M = []

		for i in range(generator.nrows):
			M.append(g_mul_perm.T.solve(generator[i])[1])

		return matrix.from_vectors(M)


	def attack(self, pub_key):

		is_dual_code = self.m <= 2 * self.r

		g_for_rm_d = pub_key

		if is_dual_code:
			self.r = self.m - 1 - self.r 
			g_for_rm_d = g_for_rm_d.orthogonal


		(a, b, d) = self.extended_gcd(self.m - 1, self.r)
		
		rm_d = self.generate_rm_d(g_for_rm_d, a, b, d)

		rm_1 = rm_d

		if d != 1:

			rm_d_minus_1 = ms.MinderShokrollahi(d, self.m).attack(rm_d)

			rm_1 = self.circle_dot_prod(rm_d.orthogonal, rm_d_minus_1).orthogonal

		P = self.find_permutation(rm_1)

		if is_dual_code:
			self.r = self.m - 1 - self.r

		M = self.find_nonsingular(pub_key, rm.generator(self.r, self.m) * P)

		return (M, P)

	def check(self, pub_key, M, P):
		return M * rm.generator(self.r, self.m) * P == pub_key 


r = int(input("Enter r: "))
m = int(input("Enter m: "))


obj = ChizhovBorodin(r, m)
pub_key = obj.pubkey_gen()
print ("Pubkey:", pub_key)

(M, P) = obj.attack(pub_key)

if obj.check(pub_key, M, P):
	print ("Success!")

else: 
	print ("Fail!")

