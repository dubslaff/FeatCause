import argparse
import logging
import os
import sys
import time
import random
import pandas as pd
from pathlib import Path

from pyeda.inter import *
from pyeda.boolalg.expr import Complement, AndOp

import fm_parser
import utils
import causality

class bc:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKCYAN = '\033[96m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'
	DARK = '\033[2;37m'

def print_num(n, color=bc.ENDC):
	s = ""
	if n == 0:
		s += bc.DARK
	else:
		s += color
	s += f"{n:6.4f}" + bc.ENDC
	return s

def print_dnf(clist, kind, idict1=dict(), idict2=dict()):
	with open(tmp_path + "/" + kind.replace(" ", "_") + "." + str(i), "w") as f:
		if len(clist) > 0:
			# print the causes
			logger.info(kind + " computed (" + bc.OKGREEN + str(len(clist)) + bc.ENDC + ")")
			for e_c in clist:
				if print_flag:
					if e_c.is_one():
						logger.info("\t" + bc.WARNING + "No results (1)" + bc.ENDC)
					else:
						ostr = "\t"
						if e_c in idict1.keys():
							ostr += "b: " + print_num(idict1[e_c]) + " "
						if e_c in idict2.keys():
							ostr += "w: " + print_num(idict2[e_c]) + " "
						if e_c in idict1.keys() and e_c in idict2.keys():
							ostr += "w*b: " + print_num(idict2[e_c]*idict1[e_c]) + " "
						ostr += bc.OKGREEN + str(e_c) + bc.ENDC
						logger.info(ostr)
				f.write(str(e_c)+"\n")
		else:
			logger.info("\t" + bc.WARNING + "No results (0)" + bc.ENDC)

# define input arguments for FeatCause
arg_parser = argparse.ArgumentParser(
	description="FeatCause - compute feature causes ala ICSE'22")

base_group = arg_parser.add_argument_group("basic parameters")
base_group.add_argument("experiment", help="base experiment name ([experiment].fm for feature model, [experiment].fs for feature space, [experiment].i for effects")
base_group.add_argument("-s", "--single", default="-1", help="number of the single effect experiment to be considered")
base_group.add_argument("-fmt", "--fmtype", default="DNF", help="type of the feature model input [NO, DNF, CNF] (default: DNF)")
base_group.add_argument("-c", "--closevalid", action="store_true", help="determine feature model by (effects U non-effects)")
base_group.add_argument("-n", "--negate", action="store_true", help="switch role of effects/non-effects")
base_group.add_argument("-na", "--noatomic", action="store_true", help="do not compute atomic causes")
base_group.add_argument("-inc", "--incnum", default="0", help="incremental adding valid configurations to effects/non-effects [number of steps, 0=no]")

expl_group = arg_parser.add_argument_group("explication parameters")
expl_group.add_argument("-d", "--distribute", action="store_true", help="minimize with distribution law")
expl_group.add_argument("-m", "--minimize", action="store_true", help="compute most general cause covers")
expl_group.add_argument("-e", "--espresso", default="no", help="compute cause-effect covers with Espresso heuristics [fast, standard, signature, qm]")
expl_group.add_argument("-b", "--blame", default="no", help="compute blames [pfblame, fblame, cblame, all]")
expl_group.add_argument("-u", "--uncertainty", action="store_true", help="include uncertainty")
expl_group.add_argument("-bdd", "--bdd", action="store_true", help="return results also in BDD representation")
expl_group.add_argument("-w", "--weight", action="store_true", help="compute cause weights")
expl_group.add_argument("-i", "--interactions", action="store_true", help="compute interaction witnesses")

sanity_group = arg_parser.add_argument_group('additional sanity parameters')
sanity_group.add_argument("-sani", "--sanity", action="store_true", help="perform sanity checks")
sanity_group.add_argument("-pye", "--pyespresso", action="store_true", help="minimize with PyEDA Espresso (slow!)")
sanity_group.add_argument("-pyp", "--pyprimes", action="store_true", help="use old python way to compute primes")

info_group = arg_parser.add_argument_group('info and statistics parameters')
info_group.add_argument("-v", "--verbose", action="store_true", help="print out verbose information")
info_group.add_argument("-p", "--print", action="store_true", help="print feature causes results")
expl_group.add_argument("-pdf", "--pdf", action="store_true", help="translate results (BDDs) into PDF")
info_group.add_argument("-stat", "--statistics", action="store_true", help="write statistics to statistics file; -m and -d needed.")
info_group.add_argument("-tstat", "--timestatistics", action="store_true", help="write time statistics to statistics file; -d not needed")
args = arg_parser.parse_args()

experiment = args.experiment
enumber = args.single
fm_type = args.fmtype
blame_type = args.blame
weight_flag = args.weight
minimize_flag = args.minimize
inc_number = int(args.incnum)
bdd_flag = args.bdd
pdf_flag = args.pdf
interaction_flag = args.interactions
noatomic_flag = args.noatomic & ~(minimize_flag | weight_flag | interaction_flag)
if blame_type != "no" and noatomic_flag:
	noatomic_flag = True
negate_flag = args.negate
distribute_flag = args.distribute
noneffects_flag = args.closevalid
pyprimes_flag = args.pyprimes
sanity_flag = args.sanity
verbose_flag = args.verbose | sanity_flag
print_flag = args.print | verbose_flag  # verbose flag imposes print flag
pyespresso_flag = args.pyespresso
espresso_type = args.espresso

# initialize logger
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

# paths to the expermiment and statistic gathering results
tmp_path = "tmp/" + str(os.getpid()) + "_" + experiment.split('/')[-1]
statistic_prefix_path = "tmp/statistic_files"
if not os.path.exists(tmp_path):
	os.makedirs(tmp_path)
if not os.path.exists(statistic_prefix_path):
	os.makedirs(statistic_prefix_path)

if args.statistics:
	# initialize statistic dictionary part 1
	global statistics
	stat_file = statistic_prefix_path + '/' + experiment.split('/')[-1] + '.stats'
	if Path(stat_file).exists():
		"""
		  The statistic file related to this experiment already exists we do not need to overwrite everything
		  We overwrite only statistical data, that is calculated/outputed from our algorithm
		"""
		new_stats = False
		statistics = pd.read_csv(stat_file)
		statistics = statistics.to_dict(orient='list')
		statistics['NumberFeatures'] = []
		statistics['NumberAtomics'] = []
		statistics['NumberMinimals'] = []
		statistics['n'] = []
		statistics['k'] = []
	else:
		"""
		  There is no statistic file for this experiment, yet.
		  We have to also include statistical data to our experiment that is not produced by our algorthim
		  E.g., name of the experiment, number of valid configurations
		"""
		new_stats = True
		statistics = {}
		statistics['ExperimentID'] = []
		statistics['ExperimentName'] = []
		statistics['NumberValidConfigurations'] = []
		statistics['NumberFeatures'] = []
		statistics['EffectSetSize'] = []
		statistics['NumberAtomics'] = []
		statistics['NumberMinimals'] = []
		statistics['n'] = []
		statistics['k'] = []
elif args.timestatistics:
	if negate_flag:
		stat_file= statistic_prefix_path + '/' + experiment.split('/')[-1] + 'NegatedStats.csv'
	else:
		stat_file= statistic_prefix_path + '/' + experiment.split('/')[-1] + 'Stats.csv'
	if Path(stat_file).exists():
		statistics = pd.read_csv(stat_file)
		statistics = statistics.to_dict(orient='list')
		statistics['TimeAtomicMinimals'] = []
	else:
		raise ValueError('First run with statistics and -d')


fileFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
fileHandler = logging.FileHandler(tmp_path + "/debug.log", mode="w")
fileHandler.setFormatter(fileFormatter)
logger.addHandler(fileHandler)

consoleFormatter = logging.Formatter("%(message)s")
consoleHandler = logging.StreamHandler(sys.stdout)
consoleHandler.setFormatter(consoleFormatter)
logger.addHandler(consoleHandler)

# print some welcome stuff
logger.info(bc.HEADER + "--- FeatCause ---" + bc.ENDC)

# get starting time
total_start_time = time.time()

# the actual program
filename_univ = experiment+".fs"
feats = fm_parser.read_fs_file(filename_univ)

b_valid_feats = expr2bdd(expr(True))
univ = set()
filename_valid_feats = experiment+".fm"
if fm_type != "NO" and not os.path.isfile(filename_valid_feats):
	logger.warning("Feature model " + filename_valid_feats + " could not be read -- continue without feature model")
else:
	if fm_type == "DNF":
		b_valid_feats, univ = fm_parser.gen_bdd_dnf_file(filename_valid_feats, feats)
	elif fm_type == "CNF":
		b_valid_feats, univ = fm_parser.gen_bdd_cnf_file(filename_valid_feats, feats)

if sanity_flag:
	if len(univ - set(feats.values())) != 0:
		logger.info("Features in variability model but not feature space " + str(univ - set(feats.values())))
	if len(set(feats.values()) - univ) != 0:
		logger.info("Features in feature space but not feature space" + str(set(feats.values())-univ))
logger.info("Feature model of " + filename_valid_feats + ": " + str(utils.get_sat_num(b_valid_feats, feats)) + " configurations")

logger.info("\tPID: "+str(os.getpid()))
logger.info("\tNumber of Features: "+str(len(feats)))
logger.info(f"--- Time to read feature model: {time.time() - total_start_time:10.4f}s")

i = 0
start_iteration = True
if enumber != "-1":
	i = int(enumber)
while os.path.isfile(experiment + "." + str(i)):
	# read effects
	filename_effects = experiment + "." + str(i)
	logger.info(bc.BOLD + bc.UNDERLINE + "--- start iteration " + str(i) + " ---" + bc.ENDC)

	start_time = time.time()
	b_effects, univ = fm_parser.gen_bdd_dnf_file(filename_effects, feats)
	if negate_flag:
		b_effects = ~b_effects
	if not noneffects_flag:
		b_effects = b_effects & b_valid_feats
	n_effects = utils.get_sat_num(b_effects, feats)
	if verbose_flag:
		logger.info(f"--- Time to build effect BDD: {time.time() - start_time:10.4f}s")
		logger.info(str(n_effects) + " effects")
		if sanity_flag:
			logger.info(f"effect feats-univ: {set(feats.values()) - univ}, univ-feats: {univ - set(feats.values())}")

	# read non-effects
	start_time = time.time()
	filename_neffects = experiment + ".n" + str(i)
	univ
	if noneffects_flag:
		b_neffects, univ = fm_parser.gen_bdd_dnf_file(filename_neffects, feats)
		b_neffects = b_neffects & b_valid_feats
	else:
		b_neffects = b_valid_feats & ~b_effects
	if verbose_flag:
		logger.info("--- Time to build BDD for valid non-effects: " + "{:10.4f}".format(time.time() - start_time) + " seconds")
	n_neffects = utils.get_sat_num(b_neffects, feats)
	if verbose_flag:
		logger.info(str(n_neffects) + " valid non-effects")
		if sanity_flag:
			logger.info("valid non-effects ... feats-univ: " + str(set(feats.values()) - univ) + " univ-feats: " + str(univ - set(feats.values())))
			if not noneffects_flag:
				logger.info("sanity check: " + str(n_effects) + " valid effects")
			logger.info("sanity check (b_ne & b_e is zero): " + str((b_neffects & b_effects).is_zero()))
	if args.statistics:
		statistics['NumberFeatures'].append(len(feats))
	if args.statistics and new_stats:
		statistics['EffectSetSize'].append(n_effects)
		statistics['NumberValidConfigurations'].append(n_effects + n_neffects)
		statistics['ExperimentID'].append(i)
		statistics['ExperimentName'].append(experiment.split('/')[-1])

	if args.uncertainty:
		# read uncertainty bdd
		filename_ueffects = experiment + "UncertainConfigurations." + str(i)
		b_ueffects, _ = fm_parser.gen_bdd_dnf_file(filename_ueffects, feats)
		b_effects &= ~b_ueffects
		b_neffects &= ~b_ueffects

	# variant with primes
	iteration_start_time = time.time()
	if not noatomic_flag:
		logger.info("Compute causes via prime implicants...")

	E_causes = set()
	E_mcauses = set()

	### stepwise adding of knowledge about effect/non-effect in valid feats (for precauses)
	if inc_number > 0: # compute stepwise causes from shuffled valid feats
		if n_effects > 0 and n_neffects > 0: # if there can be causes at all
			E_effects = list(map(utils.sat_dict_to_expr, utils.get_all_sat(b_effects, set(feats.values()))))
			E_neffects = list(map(utils.sat_dict_to_expr, utils.get_all_sat(b_neffects, set(feats.values()))))
			fscores = [0] * (n_effects+n_neffects)
			precisions = [0] * (n_effects+n_neffects)
			print("Start incremental computation for", i, "with", n_effects, n_neffects)
			for j in range(inc_number):
				print(j, end=".")
				# generate shuffled list of effects and initial BDD
				random.shuffle(E_effects)
				ei = 0
				bi_effects = expr2bdd(expr(0)) # set initial BDD to false
				# generate shuffled list of non-effects and initial BDD
				random.shuffle(E_neffects)
				nei = 1
				bi_neffects = expr2bdd(E_neffects[0]) # pick the first non-effect (to have at least one counterfactual)
				# with open(f"{tmp_path}/inclist_{j}.{i}", "w") as f: # for reproducibility: add order
				# 	f.write("\n".join(map(str, E_effects)) + "\n" + "\n".join(map(str, E_neffects)))
				
				# start incremental session
				while ei+nei < n_effects+n_neffects:
					# add new sample
					if ei/nei < n_effects/n_neffects and ei < n_effects: # ratio less for effects, add effect
						bi_effects |= expr2bdd(E_effects[ei])
						ei += 1
					else: # otherwise add counterfactual
						bi_neffects |= expr2bdd(E_neffects[nei])
						nei += 1

					# E_causes, E_mcauses = causality.get_ternary_causes(bi_effects, bi_neffects, b_valid_feats, feats, tmp_path, expr=pyprimes_flag, minimization=minimize_flag)
					E_causes, E_mcauses = causality.get_causes(bi_effects, bi_neffects, feats, tmp_path, expr=pyprimes_flag, minimization=minimize_flag)
					# print([list(x)[0] for x in E_mcauses])
					# compute most certain covering -- sort by ratio uncertain/covered effects
					b_causes = expr2bdd(expr(0))
					mprecauses = list()
					while not (bi_effects & ~b_causes).is_zero():
						clist = list()
						for e_c in [list(x)[0] for x in E_mcauses]:
							b_c = expr2bdd(e_c)
							b_add = b_c & bi_effects & ~b_causes
							if not b_add.is_zero(): # add e_c as candidate to cover more incremental effects
								# covered_effects = utils.get_sat_num(b_c & (bi_effects | b_causes | ~b_valid_feats), feats)
								# hard number of additionally covered effects -- utility
								utility_ieffects = utils.get_sat_num(b_add, feats)
								# additional effects to be assumed -- costs
								cost_ieffects = utils.get_sat_num(b_c & (b_valid_feats & ~bi_effects & ~b_causes), feats)
								# covered configs with already added precauses as effects
								# covered_configs = utils.get_sat_num(b_c & (bi_effects | b_causes | ~b_valid_feats), feats)
								if isinstance(e_c, AndOp) and e_c.depth == 1:
									small = len(e_c.xs)
								else:
									small = 1
								# covered possible effects
								# covered_peffects = utils.get_sat_num(b_c & b_valid_feats, feats) /(covered_peffects+0.0001)
								clist.append((e_c,utility_ieffects/(cost_ieffects+0.0001),small))
						clist.sort(key=lambda a: a[2],reverse=False) # first sort on most covered effects
						clist.sort(key=lambda a: a[1],reverse=True) # then sort on most certain
						
						mprecauses.append(clist[0][0]) # add first cause expressions
						b_causes |= expr2bdd(clist[0][0]) # add covered effects to cover
					# print(mprecauses)

					# fscore = 2TP / (2TP + FP + FN)
					tp = utils.get_sat_num(b_causes & b_effects, feats)
					fp = utils.get_sat_num(b_causes & b_neffects, feats)
					fn = utils.get_sat_num(~b_causes & b_effects, feats)
					# print(ei, nei, 'fscore:', str(2*tp/(2*tp+fp+fn))[:5], 'precision:', str(tp/(tp+fp))[:5])
					fscores[ei+nei-1] += (2*tp/(2*tp+fp+fn))/inc_number
					precisions[ei+nei-1] += (tp/(tp+fp))/inc_number
			# write out the average fscore
			with open(f"{tmp_path}/fscore.{i}", "w") as ff:
				for j in range(n_effects+n_neffects):
					ff.write(f"{j}\t{fscores[j]}\n")
			# write out the average precision
			with open(f"{tmp_path}/precision.{i}", "w") as fp:
				for j in range(n_effects+n_neffects):
					fp.write(f"{j}\t{precisions[j]}\n")
		
	elif not noatomic_flag:
		E_causes, E_mcauses = causality.get_causes(b_effects, b_neffects, feats, tmp_path, expr=pyprimes_flag, minimization=minimize_flag)

	if bdd_flag:
		b_causes = expr2bdd(expr(0))
		for e_c in E_causes:
			b_causes |= expr2bdd(e_c)
		with open(f"{tmp_path}/cause_bdd{i}.dot", "w") as f:
			f.write(b_causes.to_dot())
		if pdf_flag:
			os.system(f"dot -Tpdf {tmp_path}/cause_bdd{i}.dot > {tmp_path}/cause_bdd{i}.pdf")
		if minimize_flag:
			j = 0
			for E_mc in E_mcauses:
				b_mcauses = expr2bdd(expr(0))
				for e_c in E_mc:
					b_mcauses |= expr2bdd(e_c)
				with open(f"{tmp_path}/mcause_bdd{i}_{j}.dot", "w") as f:
					f.write(b_mcauses.to_dot())
				if pdf_flag:
					os.system(f"dot -Tpdf {tmp_path}/mcause_bdd{i}_{j}.dot > {tmp_path}/mcause_bdd{i}_{j}.pdf")
				j += 1

	B_causes = dict()
	if blame_type == "cblame" or blame_type == "all":
		logger.info("Compute partial blame of causes...")
		for e_c in E_causes:
			blame = causality.get_partial_blame(e_c, E_causes, b_effects, b_neffects, set(feats.values()))
			B_causes[e_c] = blame

	W_causes = dict()
	if weight_flag:
		logger.info("Compute weight of causes...")
		W_causes = causality.get_cause_weights(E_causes, b_effects, feats)

	if not noatomic_flag:
		print_dnf(E_causes, "O Atomic feature causes", B_causes, W_causes)
		if distribute_flag and len(E_causes)>0:
			cexpr = Or(*E_causes)
			dexpr = utils.distribute(cexpr)
			pstr = "distributive law simplification"
			if print_flag:
				pstr = bc.FAIL + str(dexpr) + bc.ENDC
			logger.info("X\tDLS [" + str(cexpr.size) + " -> " + str(dexpr.size) + "]: " + pstr)
			if args.statistics:
				statistics['n'].append(cexpr.size)
				statistics['k'].append(dexpr.size)
		if distribute_flag and len(E_causes) == 0:
			if args.statistics:
				statistics['n'].append(0)
				statistics['k'].append(0)
			if sanity_flag:
				logger.info("sanity check (reduced expression semantically equivalent?): " + bc.OKBLUE + str(cexpr.equivalent(dexpr)) + bc.ENDC)
		if args.statistics:
			statistics['NumberAtomics'].append(len(E_causes))

	# cross-check with minimized semantics
	# primes
	if minimize_flag:
		with open(tmp_path + "/most_general_causes."+str(i), "w") as f:
			if len(E_mcauses) > 0:
				E_all_mcauses = list()
				logger.info("X Most general causes (" + bc.OKCYAN + str(len(E_mcauses)) + bc.ENDC + ")")
				if args.statistics:
					statistics['NumberMinimals'].append(len(E_mcauses))
				### actually we need not to recompute blame, but because of different keys in mcauses/causes - TODO: more elegant approach...
				BM_causes = dict()
				if blame_type == "cblame" or blame_type == "all":
					logger.info("Compute partial blame of most general causes...")
					for E_c in E_mcauses:
						l = sorted(E_c, key=utils.get_expr_length)
						blame = causality.get_partial_blame(l[0], E_causes, b_effects & b_valid_feats, b_neffects & b_valid_feats, set(feats.values()))
						BM_causes[l[0]] = blame

				for E_c in E_mcauses:
					E_all_mcauses += E_c
					l = sorted(E_c, key=utils.get_expr_length)
					slist = map(str, l)
					s = bc.OKCYAN + " | ".join(slist) + bc.ENDC
					s = s.replace( " | ", bc.WARNING + " | " + bc.OKCYAN)
					if print_flag:
						b = ""
						if len(BM_causes) > 0:
							b = "b: " + print_num(BM_causes[l[0]]) + " "
						logger.info("X\t" + b + s)
					f.write("\n")
				if distribute_flag:
					cexpr = Or(*E_all_mcauses)
					dexpr = utils.distribute(cexpr)
					pstr = "distributive law simplification"
					if print_flag:
						pstr = bc.FAIL + str(dexpr) + bc.ENDC
					logger.info("X\tDLS [" + str(cexpr.size) + " -> " + str(dexpr.size) + "]: " + pstr)
					if sanity_flag:
						logger.info("sanity check (reduced expression semantically equivalent?): " + bc.OKBLUE + str(cexpr.equivalent(dexpr)) + bc.ENDC)
			else:
				logger.info("X\t" + bc.WARNING + "No most general causes" + bc.ENDC)
				if args.statistics:
					statistics['NumberMinimals'].append(0)
		if args.timestatistics:
			statistics['TimeAtomicMinimals'].append(f"{time.time() - iteration_start_time:10.4f}")
		logger.info(f"--- Time for atomic and most general causes: {time.time() - iteration_start_time:10.4f}s")
	elif not noatomic_flag:
		logger.info(f"--- Time for atomic causes: {time.time() - iteration_start_time:10.4f}s")

# feature interactions
	if interaction_flag:
		if len(E_causes)>0:
			t, E_iw = causality.get_tway_interactions(E_causes)
			print_dnf(E_iw, "I "+str(t)+"-way interaction witnesses")
		else:
			logger.info("I No interaction witnesses since no atomic causes.")

# blame computations
	if "blame" in blame_type or blame_type == "all":
		if len(E_causes) > 0:
			iteration_start_time = time.time()
			if blame_type == "pfblame" or blame_type == "all":
				logger.info("R Partial blame computations ...")
				blame_candidates = fm_parser.read_fr_file(experiment + ".fr")
				for candidate in blame_candidates:
					e_l = expr("And(" + str(candidate) + ")")
					partial_blame = causality.get_partial_blame(e_l, E_causes, b_effects, b_neffects, set(feats.values()))
					with open(tmp_path + "/partial_blame.csv", "a") as f:
						f.write(str(candidate) + "\t" + str(partial_blame) + "\n")
						print(f"R {candidate}: \t {print_num(partial_blame)}")
			if blame_type == "fblame" or blame_type == "all":
				logger.info("R Blame computations for each feature")
				with open(tmp_path + "/feature_blame.csv", "a") as f:
					rstr = str(i)
					istr = "# nr"
					r0 = "none"
					for x in sorted(list(feats.values()), key=str):
						r0 = causality.get_partial_blame(expr(~x), E_causes, b_effects, b_neffects, set(feats.values()))
						rstr += "\t" + str(r0)
						r1 = causality.get_partial_blame(expr(x), E_causes, b_effects, b_neffects, set(feats.values()))
						rstr += "\t" + str(r1)
						if start_iteration:
							istr += "\t" + str(~x) + "\t" + str(x)
						logger.info("R\t[" + print_num(r0) + "," + print_num(r1) +"] " + str(x))
					if start_iteration:
						f.write(istr+"\n")
						start_iteration = False
					f.write(rstr+"\n")
			logger.info(f"--- Time for blame computation: {time.time() - iteration_start_time:10.4f}s")
		else:
			logger.info("R No blame computation since no causes")

	# SOP minimization with Espresso
	if espresso_type != "no":
		iteration_start_time = time.time()
		logger.info("E Compute a (nearly) minimal SOP via Espresso...")
		e = utils.espresso_bdd(b_effects, ~b_effects & ~b_neffects, feats, tmp_path, type=espresso_type)
		print_dnf(e, "E Espresso minterms")
		logger.info(f"E --- Time for Espresso: {time.time() - iteration_start_time:10.4f}s")

	# SOP minimization with Espresso via PyEDA
	if pyespresso_flag:
		iteration_start_time = time.time()
		logger.info("Y Compute a nearly minimal SOP via Espresso...")
		msop = utils.pyespresso_bdd(b_effects, ~b_effects & ~b_neffects)
		if isinstance(msop, Complement) or isinstance(msop, Variable) or isinstance(msop, AndOp):
			logger.info("Y Espresso (" + bc.OKBLUE + "1" + bc.ENDC + ")")
			if print_flag:
				logger.info("Y\t" + bc.OKBLUE + str(msop) + bc.ENDC)
		elif msop.is_zero():
			logger.info("Y Espresso (0): " + bc.WARNING + "False" + bc.ENDC)
		elif msop.is_one():
			logger.info("Y Espresso (0): " + bc.WARNING + "True" + bc.ENDC)
		else:
			logger.info("Y Espresso (" + bc.OKBLUE + str(len(msop.xs)) + bc.ENDC + ")")
			if print_flag:
				for e_c in msop.xs:
					logger.info("Y \t" + bc.OKBLUE + str(e_c) + bc.ENDC)
		logger.info(f"Y --- Time for PyEspresso: {time.time() - iteration_start_time:10.4f}s")

	# increment iteration
	if enumber != "-1":
		break
	else:
		i += 1

if args.statistics or args.timestatistics:
	stat_outdir = statistic_prefix_path
	if not os.path.exists(stat_outdir):
		os.mkdir(stat_outdir)
	if negate_flag:
		print(statistics)
		dataframe = pd.DataFrame.from_dict(statistics)
		dataframe.to_csv(stat_outdir + '/' + experiment.split('/')[-1] + 'NegatedStats.csv', index=False)
		print(dataframe)
	else:
		print(statistics)
		dataframe = pd.DataFrame.from_dict(statistics)
		dataframe.to_csv(stat_outdir + '/' + experiment.split('/')[-1] + 'Stats.csv', index=False)
		print(dataframe)


##########

# finalize output
logger.info("------------- end iterations")
logger.info(f"Total time for computation: {time.time() - total_start_time:10.4f}s")
