from pyeda.inter import *
import utils

LO_PROD = 0  # list of products
LO_POSPROD = 1  # list of positive products

def read_fs_file(filename):
	with open(filename, "r") as f:
		univ = set()
		for line in f:
			line = line.replace("\n", "")
			# skip comments and empty lines
			if line.startswith("#") or len(line) == 0 or line == " ":
				continue
			# get the universe of variables used
			e = expr(line)
			univ |= set(e.inputs)
	# convert to a dictionary
	i = 0
	duniv = dict()
	for u in sorted(list(univ)): # sort the list for reproducibility
		duniv[i] = u
		i += 1
	return duniv

def read_dnf_file(filename, feats, itype=LO_PROD):
	with open(filename, "r") as f:
		slist = list()
		univ = set()
		for line in f:
			line = line.replace("\n", "")
			# skip comments and empty lines
			if line.startswith("#") or len(line) == 0 or line == " ":
				continue
			slist.append(line)
			# compute the universe of variables used
			e = expr(line)
			univ |= set(e.inputs)

		if len(slist) == 0:
			result = expr(0) # empty disjunction is false
		elif len(slist) == 1:
			result = expr(slist[0]) # single disjunction is expression
		else:
			# determine the expression
			if itype == LO_PROD:  # input is in list of conjunctions
				rs = "Or(" + ", ".join(slist) + ")"
				# rs = "(" + ") | (".join(slist) + ")"
			elif itype == LO_POSPROD:  # input is in a list of products, i.e., only comprising positive literals
				# go through all expressions and complete the products
				nlist = list()
				for i in range(len(slist)):
					missing = feats - set(expr(slist[i]).inputs)
					n = slist[i]
					for m in missing:
						n += ' & ~' + str(m)
					nlist.append(n)
				rs = "(" + ") | (".join(nlist) + ")"
			result = expr(rs)
		return result, univ

def gen_bdd_dnf_file(filename, feats):
	with open(filename, "r") as f:
		result = expr2bdd(expr(0))
		univ = set()
		for line in f:
			line = line.replace("\n", "")
			# skip comments and empty lines
			if line.startswith("#") or len(line) == 0 or line == " ":
				continue
			e = expr(line)
			b = expr2bdd(e)
			# print("expression: ", e)
			result |= b
			# print(utils.get_sat_num(result, feats), bdd2expr(result))
			# compute the universe of variables used
			univ |= set(e.inputs)
		return result, univ

def gen_bdd_cnf_file(filename, feats):
	with open(filename, "r") as f:
		result = expr2bdd(expr(1))
		univ = set()
		j = 0
		for line in f:
			j += 1
			print("{:5}".format(j), end="\r")
			line = line.replace("\n", "")
			# skip comments and empty lines
			if line.startswith("#") or len(line) == 0 or line == " ":
				continue
			e = expr(line)
			# print("expression: ", e)
			result &= expr2bdd(e)
			# print(utils.get_sat_num(result, feats), bdd2expr(result))
			# compute the universe of variables used
			univ |= set(e.inputs)

		return result, univ

def read_dimacs_file(filename, feats):
	with open(filename, "r") as f:
		univ = set()
		clist = list()
		k = 0
		for line in f:
			# replace each 0 at the end of the line
			line = line.replace(" 0\n", "")
			line = line.replace("\n", "")
			# skip comments and empty lines
			if line.startswith("c") or line.startswith("p") or len(line) == 0 or line == " ":
				continue
			spl = line.split(" ")
			dlist = list()
			for i in range(len(spl)):
				s = spl[i]
				if s.startswith("-"):
					s = "Not(x" + s[1:] + ")"
				else:
					s = "x" + s
				dlist.append(s)
			nline = "Or("+",".join(dlist)+")"
			print(k, nline)
			k += 1
			# compute the universe of variables used
			e = expr(nline)
			clist.append(e)
			univ |= set(e.inputs)
		result = And(*clist)
		return result, univ

def read_fr_file(filename):
	with open(filename, "r") as f:
		r_tuples = []
		for line in f:
			line = line.replace("\n", "")
			# skip comments and empty lines
			if line.startswith("#") or len(line) == 0 or line == " ":
				continue
			r_tuples.append(line)
	return r_tuples
