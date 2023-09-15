import random
import subprocess
import os
import time
import utils

from pyeda.boolalg.espresso import espresso, set_config
from pyeda.boolalg.minimization import CONFIG, _cover2exprs
from pyeda.boolalg.expr import *
from pyeda.inter import *


###############################################
### feature causality
###############################################
def get_causes(b_e, b_ne, umap, tmp_path, expr=False, minimization=True):
	""" takes BDDs {b_e} for effects and valid non-effects {b_ne} and returns a 
		set of causes as Boolean formulas using prime-implicant computation
	"""
	E_pfcauses = set()
	E_mpfcauses = set()
	B_clist = list()
	B_eb = list()
	N_eb = list()

	b = ~b_ne

	if not b.is_zero() and not b.is_one():
		# print("O\tCompute primes...")
		if expr:
			primes = utils.get_expr_primes(b)
		else:
			primes = utils.espresso_bdd(b, None, umap, tmp_path, type="primes")

		# primes.sort(key=utils.get_andsize, reverse=False)
		j = 0
		start_time = time.time()
		for e_a in primes: # go through all the prime implicants and add if there is a valid evaluation
			b_a = expr2bdd(e_a)
			b_ea = b_a & b_e
			if j > 200:
				print(f"{j:5}/{len(primes)} primes checked, {len(B_clist):5}/{len(E_pfcauses):5} mins/atomics [{time.time() - start_time:10.4f}s]", end="\r")
				if j%int(len(primes)/20+200) - 199 == 1:
					print() # make a log permanent to compare
			j += 1
			if not b_ea.is_zero():
				E_pfcauses.add(e_a) # add cause to list
				if minimization:
					n_ea = utils.get_sat_num(b_ea, umap)
					# check whether covered - go through equivalence classes
					for i in reversed(range(len(B_clist))):
						b_b = B_clist[i][0]
						b_eb = B_eb[i]
						n_eb = N_eb[i]
						if (n_ea <= n_eb) and (b_ea & ~b_b).is_zero(): # b_a covered by existing b_b
							# if len(b_b.support) < len(b_a.support): # and existing are smaller -> break
							# 	break
							# if len(b_b.support) == len(b_a.support) and (n_ea == n_eb) and ((b_eb & ~b_a).is_zero()):
							if (n_ea == n_eb) and ((b_eb & ~b_a).is_zero()):
								# b_a covers b_b -> add to equivalence class
								B_clist[i].append(b_a)
							break
						elif (n_ea >= n_eb) and ((b_eb & ~b_a).is_zero()): # if b_b is covered by b_a -> remove b_b
							B_clist.pop(i)
							B_eb.pop(i)
							N_eb.pop(i)
					else:
						B_clist.append([b_a])
						B_eb.append(b_ea)
						N_eb.append(n_ea)
	print(80*" ", end="\r") # dirty flush line
	for C in B_clist:
		c = set()
		for b in C:
			c.add(bdd2expr(b))
		E_mpfcauses.add(frozenset(c))
	return E_pfcauses, E_mpfcauses

def get_ternary_causes(b_e, b_ne, b_all, umap, tmp_path, expr=False, minimization=True):
	""" takes BDDs {b_e} for effects, non-effects {b_ne}, and valids {b_all} and returns a 
		set of causes as Boolean formulas using prime-implicant computation
	"""
	E_pfcauses = set()
	E_mpfcauses = set()
	B_clist = list()
	B_eb = list()
	N_eb = list()

	b = ~b_all | b_e

	if not b.is_zero() and not b.is_one():
		# print("O\tCompute primes...")
		if expr:
			primes = utils.get_expr_primes(b)
		else:
			primes = utils.espresso_bdd(b, b_all & ~b_e & ~b_ne, umap, tmp_path, type="primes")

		# primes.sort(key=utils.get_andsize, reverse=False)
		j = 0
		start_time = time.time()
		for e_a in primes: # go through all the prime implicants and add if there is a valid evaluation
			b_a = expr2bdd(e_a)
			b_ea = b_a & b_e
			if j > 200:
				print(f"{j:5}/{len(primes)} primes checked, {len(B_clist):5}/{len(E_pfcauses):5} mins/atomics [{time.time() - start_time:10.4f}s]", end="\r")
				if j%int(len(primes)/20+200) - 199 == 1:
					print() # make a log permanent to compare
			j += 1
			if not b_ea.is_zero():
				E_pfcauses.add(e_a) # add cause to list
				if minimization:
					n_ea = utils.get_sat_num(b_ea, umap)
					# check whether covered - go through equivalence classes
					for i in reversed(range(len(B_clist))):
						b_b = B_clist[i][0]
						b_eb = B_eb[i]
						n_eb = N_eb[i]
						if (n_ea <= n_eb) and (b_ea & ~b_b).is_zero(): # b_a covered by existing b_b
							# if len(b_b.support) < len(b_a.support): # and existing are smaller -> break
							# 	break
							# if len(b_b.support) == len(b_a.support) and (n_ea == n_eb) and ((b_eb & ~b_a).is_zero()):
							if (n_ea == n_eb) and ((b_eb & ~b_a).is_zero()):
								# b_a covers b_b -> add to equivalence class
								B_clist[i].append(b_a)
							break
						elif (n_ea >= n_eb) and ((b_eb & ~b_a).is_zero()): # if b_b is covered by b_a -> remove b_b
							B_clist.pop(i)
							B_eb.pop(i)
							N_eb.pop(i)
					else:
						B_clist.append([b_a])
						B_eb.append(b_ea)
						N_eb.append(n_ea)
	print(80*" ", end="\r") # dirty flush line
	for C in B_clist:
		c = set()
		for b in C:
			c.add(bdd2expr(b))
		E_mpfcauses.add(frozenset(c))
	return E_pfcauses, E_mpfcauses

#####################################
### Responsibility/Blame
#####################################

def get_distance(d_e, d_ne, m=1000):
	"""	Hamming distance of two Boolean vectors
		{d_e}, {d_ne}:	the two vectors as dicts
		{m}:			break-limit (to speed-up minimum computation)
	"""
	d = 0
	for y in set(d_ne.keys()) & set(d_e.keys()):
		if d_e[y] != d_ne[y]:
			d += 1
		if d >= m:
			break
	return d

def get_single_feature_responsibility(x, E_c, d_e, D_ne):
	"""	Get single responsibility of {x}
		{x}:	single feature
		{E_c}:	causes (list of expressions)
		{d_e}:	context (interpretation dictionary)
		{D_ne}:	non-effects (list of interpretation dictionaries)
	"""
	# set the minimal distance to impossible high value
	mindist = [len(d_e.keys())+1, len(d_e.keys())+1]
	for e_c in E_c:
		D_c = e_c.satisfy_all()
		d_c = next(D_c, None) # there should be only one object
		if next(D_c, None) != None:
			print("Oops, more than one element in total feature configuration")
		if x not in set(d_c.keys()):
			continue
		for y in d_c.keys():
			if d_c[y] != d_e[y]:
				break
		else: # cause has same polarity as effect as d_e[x] == d_c[x]
			for d_ne in D_ne:
				if x not in set(d_ne.keys()) or d_ne[x] != d_e[x]: # how many variables else have to switch?
					mindist[d_e[x]] = get_distance(d_e, d_ne,mindist[d_e[x]])
					# this version is from the teams paper: fix all variables in the cause, add length
					# mindist[d_e[x]] = min(mindist[d_e[x]], len(d_c) + get_distance(d_e, d_ne, d_c))
	r = [0,0]
	for i in [0,1]:
		if mindist[i] != len(d_e.keys())+1:
			r[i] = 1/mindist[i]
	return r # positive and negative polarity

def get_partial_blame(e_p, E_c, b_e, b_ne, univ):
	"""	Compute the partial interpretation blame
		{e_p}:	partial interpretation [as And(*Literal)]
		{E_c}:	a list of causes (as expressions)
		{b_e}:	BDD for effects
		{b_ne}:	BDD for non-effect
		{univ}:	set of features
	"""
	start_time = time.time()

	# get literals of e_p
	s_p = e_p.support
	if isinstance(e_p, AndOp) and e_p.depth == 1:
		l_p = e_p.xs
	elif isinstance(e_p, Literal):
		l_p = [e_p]
	else:
		print("Oops, expected partial interpretation for responsibility!")
		return 0

	# check whether there is a cause for the partial interpretation
	for e_c in E_c:
		# get literals of e_p
		if isinstance(e_c, AndOp) and e_c.depth == 1:
			l_c = e_c.xs
		elif isinstance(e_c, Literal):
			l_c = [e_c]
		else:
			print("Oops, expected partial interpretation for causes!")
			continue
		# check whether polarities in partial interpretation are the same as in the cause
		if len(set(l_p) & set(l_c)) != 0 and len(set(l_p)-set(l_c)) == 0:
			break # cause found -- compute distances
	else:
		return 0 # no cause found - responsibility 0 for all effects

	D_ne = list() # list of all total negative effects
	for dbdd in (b_ne & expr2bdd(~e_p)).satisfy_all(): # only those negative effects that could be results of switching
		d = dict()
		for x in dbdd.keys():
			d[bdd2expr(x)] = dbdd[x]
		D_ne.append(d)

	## generate a list of switched negative effects w.r.t. a single feature
	## filter out those that have not switched polarity
	D_xne = dict()
	for x in s_p: # for all variables in the support of partial interpretation
		D_xne[x] = list()
		for d_ne in D_ne:
			d_xne = dict(d_ne)
			if x in l_p: # switch the polarity of variable x - covers the case where x is not defined in d_ne
				d_xne[x] = 0
			else:
				d_xne[x] = 1
			if x in d_ne.keys() and d_ne[x] != d_xne[x]: # if the polarity from the switch is different
				continue
			D_xne[x].append(d_xne)

	D_pe = utils.get_all_sat(b_e & expr2bdd(e_p), univ) # list of dicts of all compatible total effects
	dtime = time.time()

	# minimize distances
	sum = 0
	k = 0
	for d_pe in D_pe: # d_pe is an effect that satisfies partial interpretation
		mindist = len(univ)+1 # impossible high value
		for x in s_p: # for all variables in the support of partial interpretation
			skeys = set(d_pe.keys()) - {x}
			for d_xne in D_xne[x]:
				d = 1
				for y in set(d_xne.keys()) & skeys:
					if d_pe[y] != d_xne[y]:
						d += 1
						if d >= mindist:
							break
				else:
					mindist = d
				if mindist == 1:
					break # the smallest value possible
			if mindist == 1:
				break # the smallest value possible
		k += 1
		print("\t{:5}".format(k), "/", len(D_pe), end="\r")
			
		if mindist != len(univ)+1:
			sum += 1/mindist
		
	u = utils.get_sat_num(b_e, univ)
	r = sum/u
	return r

#####################################
### Cause weights
#####################################

def get_cause_weights(E_c, b_e, umap):
	# reverse name to position for pla string
	imap = {expr(v): k for k, v in umap.items()}

	D_e = utils.get_all_sat(b_e, set(umap.values())) # list of dicts of all compatible total effects
	# print(D_e)
	P_e = list()
	for d_e in D_e:
		# convert to string
		p_e = utils.clause2pla(utils.sat_dict_to_expr(d_e), imap)[:-2]
		P_e.append(p_e)

	# go through all causes and count number of covers
	Q_c = dict()
	for p_e in P_e:
		Q_c[p_e] = 0
		for e_c in E_c:
			p_c = utils.clause2pla(e_c, imap)[:-2]
			if utils.is_covered(p_c,p_e):
				Q_c[p_e] += 1
	# print(Q_c)

	# go through all causes and build up fractions
	W_c = dict()
	for e_c in E_c:
		W_c[e_c] = 0
		for p_e in P_e:
			p_c = utils.clause2pla(e_c, imap)[:-2]
			if utils.is_covered(p_c,p_e):
				W_c[e_c] += 1/Q_c[p_e]
	# print(W_c)

	return W_c

#####################################
### Feature Interactions
#####################################

def get_tway_interactions(E_c):
	if len(E_c) > 0:
		t = min([utils.get_expr_length(c) for c in E_c])-1
		E_iw = {c for c in E_c if utils.get_expr_length(c)-1 == t}
		t = max(1, t)
		return t, E_iw
	else:
		return 0, E_c