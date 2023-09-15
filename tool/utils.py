import random
import subprocess
import os
import time

from pyeda.boolalg.espresso import espresso, set_config
from pyeda.boolalg.minimization import CONFIG, _cover2exprs
from pyeda.boolalg.expr import *
from pyeda.inter import *

##############
# SAT, EXPR
##############

def get_sat_num(b, univ):
	""" takes a BDD {b}, a set of variables {univ}
		returns the number of satisfying evaluations
	"""
	n = 0
	for x in b.satisfy_all():
		n += pow(2, len(univ)-len(x))
	return n

def complete_sat_dict(d, univ):
	""" takes a dict with 0/1 assignments and completes it with missing possibilities
		with respect to a set (!) of variables
		returns a set of such completions
	"""
	rset = list()
	euniv = list(univ - d.keys())
	if len(euniv) == 0:
		rset.append(d)
	else:
		nd = d.copy()
		x = euniv[0]
		nd[x] = 0
		rset.extend(complete_sat_dict(nd, univ))
		nd = d.copy()
		nd[x] = 1
		rset.extend(complete_sat_dict(nd, univ))
	return rset

def sat_dict_to_expr(d):
	"""	translates a 0/1 dictionary {d} to an equivalent conjunction
	"""
	elist = list()
	for x in d.keys():
		if d[x] == 1:
			elist.append(x)
		else:
			elist.append(~x)
	return And(*elist)

def get_all_sat(b, univ):
	"""	get all satisfying assignments of a BDD {b}
	"""
	slist = list()
	for dbdd in b.satisfy_all():
		d = dict()
		# need to map BDD variables and EXPR variables
		for x in dbdd.keys():
			d[bdd2expr(x)] = dbdd[x]
		slist.extend(complete_sat_dict(d, univ))
	return slist

##############
# Distributive law simplification
##############
def get_expr_length(e):
	"""	simply return length of an expression {e} (e.g. for custom sorting)
	"""
	return e.size

def distribute(ex, exhaustive=False, s_sorted=False):
	"""	takes an expression {ex} and performs DLS heuristics
		{exhaustive}:	if set to true, performs an exhaustive search for
						minimal length DLS factorization (slow!)
		{s_sorted}:		if set to true, factor out lengthy expressions first
	"""
	cand = ex
	if isinstance(ex, OrOp) and len(ex.xs) > 1:
		if s_sorted:
			sex = sorted(ex.xs, key=get_expr_length)
		else:
			sex = ex.xs
		histo = dict()
		# build histogram of contained variables
		for x in sex:
			if isinstance(x, AndOp) and len(x.xs) > 1:
				for y in x.xs:
					if y in histo.keys():
						histo[y] += 1
					else:
						histo[y] = 1
		# factorize the most frequent variable
		s = sorted(histo, key=histo.get, reverse=True)
		if len(s) == 0:
			return cand
		if not exhaustive:
			s = list(filter(lambda x: x == s[0], s))

		for alpha in s:
			p0_list = list()
			p1_list = list()
			for x in sex:
				if isinstance(x, AndOp) and len(x.xs) > 1:
					if alpha in x.xs:
						newxs = list(x.xs)
						newxs.remove(alpha)
						p1_list.append(And(*newxs))
						continue
				p0_list.append(x)
			if len(p0_list) == 0:
				newcand = And(alpha,distribute(Or(*p1_list)))
			else:
				newcand = Or(distribute(Or(*p0_list)), And(alpha,distribute(Or(*p1_list))))
			if newcand.size < cand.size:
				cand = newcand
	return cand

############################
# I/O and conversion
############################
def is_covered(p_check, p_total):
	"""	check whether {p_total} is covered by {p_check}
		given as PLA-strings
	"""
	for i in range(len(p_total)):
		if p_check[i] == "-" or p_total[i] == p_check[i]:
			continue
		else:
			return False
	return True

def pla_bdd(outf, b_fun, b_dc, umap):
	""" PLA output to run external espresso with
		{outf}:		PLA file to write to
		{b_fun}:	BDD of the ON-function
		{b_dc}:		BDD of the don't cares (DCs)
		{umap}:		dict from variable integers to variables (variable universe)
	"""
	pla_fun = bdd2pla(b_fun, umap, ftype="ON")
	pla_dc = bdd2pla(b_dc, umap, ftype="DC")

	with open(outf, "w") as f:
		f.write(".ftype fd\n")
		f.write(".i "+str(len(umap))+"\n")
		f.write(".o 1\n")
		f.write(pla_fun + pla_dc)


def clause2pla(e, imap):
	""" generate pla line for a clause - assumes And(*xs) input
		{e}:	And(*xs) expression
		{imap}:	dict mapping variables to its integers
	"""
	pmap = dict()

	for ex in e.iter_dfs():
		if isinstance(ex, Complement) and (~ex) in imap.keys():
			pmap[imap[~ex]] = "0"
		elif isinstance(ex, Variable) and (ex) in imap.keys():
			pmap[imap[ex]] = "1"

	pla = ""
	for i in range(len(imap)):
		if i in pmap.keys():
			pla += pmap[i]
		else:
			pla += "-"
	pla += " 1"
	return pla

def bdd2pla(b_fun, umap, ftype="ON"):
	""" function to convert BDDs to PLA string following a universe-map order
		{b_fun}:	BDD of the set
		{umap}:		mapping variable integers to variables (variable universe)
		{ftype}:	type of set described by {b_fun}, i.e., ON, OFF, or DC
	"""
	# first cover the case of not providing input (e.g. to cover empty DC)
	if b_fun == None:
		return ""
	print("Convert BDD to PLA", end="\r")
	start_time = time.time()
	pla_string = ""
	# reverse name to position for pla
	imap = {expr(v): k for k, v in umap.items()}

	if b_fun.is_one():
		pla_string = len(umap)*"-" + " 1"
	elif not b_fun.is_zero():
		e_dnf = bdd2expr(b_fun)
		print(f"DNF generated [{time.time() - start_time:10.4f}s ]", end="\r")
		clist = list()
		if isinstance(e_dnf, Literal) or isinstance(e_dnf, AndOp) or e_dnf.is_one():
			clist.append(e_dnf)
		elif isinstance(e_dnf, OrOp):
			clist = list(e_dnf.xs)

		j = 0
		for e_x in clist:
			s = clause2pla(e_x, imap)
			pla_string += s + "\n"
			j += 1
			print(f"{j:5}/{len(clist)} covers generated [{time.time() - start_time:10.4f}s ]", end="\r")

	print(60*" ", end="\r") # dirty flush line
	# adapt output depending on whether function BDD
	if ftype == "OFF":
		pla_string = pla_string.replace(" 1", " 0")
	elif ftype == "DC":
		pla_string = pla_string.replace(" 1", " -")
	return pla_string


def espresso_bdd(b_fun, b_dc, umap, tmp_path, type="standard"):
	""" Call external espresso for QM/Signature/Fast etc.
		{b_fun}:	BDD of the function
		{b_dc}:		BDD for don't cares (DC)
		{umap}:		mapping variable integers to variables (variable universe)
		{tmp_path}:	path to output intermediate results to
		{type}:		Espresso algorithm (qm,signature,fast,standard)
					or compute primes (for feature causes)
	"""
	start_time = time.time()
	inf = tmp_path + "/in.pla"
	esf = tmp_path + "/espresso.pla"
	# generate PLA from BDD
	pla_bdd(inf, b_fun, b_dc, umap)

	# print(f"O\tTime for BDD PLA export: {time.time() - start_time:10.4f}s")
	start_time = time.time()
	with open(esf, "w") as fout:
		if type == "qm":
			subprocess.Popen(["./espresso", "-Dqm", "-t", inf], stdout=fout).wait()
		elif type == "signature":
			subprocess.Popen(["./espresso", "-Dsignature", "-t", inf], stdout=fout).wait()
		elif type == "fast":
			subprocess.Popen(["./espresso", "-efast", "-Dsingle_output", "-t", inf], stdout=fout).wait()
		elif type == "standard":
			subprocess.Popen(["./espresso", "-Dsingle_output", "-t", inf], stdout=fout).wait()
		else:
			subprocess.Popen(["./espresso", "-Dprimes", "-t", inf], stdout=fout).wait()
	es_list = espr2exprlist(esf, umap)
	# print(f"O\tTime for Espresso {type}: {time.time() - start_time:10.4f}s")
	return es_list

def espr2exprlist(fname, umap):
	""" transform espresso output to an expression list
		{fname}:	file to read
		{umap}:		mapping variable integers to variables (variable universe)
	"""
	plist = list()
	with open(fname, "r") as fin:
		for line in fin:
			if line.startswith("#") or line.startswith(".") or len(line) == 0 or line == " ":
				continue
			clist = list()
			for i in range(len(umap)):
				if line[i] == "1":
					clist.append(expr(umap[i]))
				elif line[i] == "0":
					clist.append(~expr(umap[i]))
			plist.append(And(*clist))
	return plist

def get_expr_primes(b):
	""" Compute primes by PyEDA (very slow, for sanity)
		{b}:	input BDD
	"""
	plist = list()
	start_time = time.time()
	print("O\tCompute primes by PyEDA ...")
	e_impl = bdd2expr(b).complete_sum()

	if isinstance(e_impl, Literal) or e_impl.is_one():
		plist.append(e_impl)
	else:
		plist = list(e_impl.xs)

	print(f"O\tTime for {len(plist)} primes: {time.time() - start_time:10.4f}s")
	return plist
