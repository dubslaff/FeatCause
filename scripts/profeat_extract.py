# a simple hack to generate DNF effects from Prism/ProFeat logs that FeatCause understands

import argparse
import re
import shlex

vars_rex = re.compile('Variables:.*')
prop_rex = re.compile('Model checking: filter\(printall, .*')
result_rex = re.compile('[0-9]*:\(.*')
end_rex = re.compile('Time for model checking:.*')


# main program
if __name__ == '__main__':
    # define input arguments for the extractor script
    argparser = argparse.ArgumentParser(description="Prism Extract - generate DNF effects for FeatCause")
    argparser.add_argument("experiment", help="path to the experiment containing the log (.log) and feature space (.fs), e.g., ./advanced/BSN")
    args = argparser.parse_args()
    experiment = args.experiment

    # first read feature variables
    feats = set()
    with open(experiment + ".fs", "r") as l:
        for line in l:
            line = line.replace("\n", "")
            feats.add(line)

    # then the log
    with open(experiment+".log", "r") as log:
        d_feats = dict() # dictionary to map positions to features
        i = -1 # keep track of the property, first is assumed to be the feature model
        s = "" # initialize empty string
        for lline in log:
            if vars_rex.match(lline):
                # build directory of feature variables
                lline = lline.replace(" \n", "")
                lline = lline.replace("Variables:   ", "")
                v = lline.split(" ")
                for k in range(len(v)):
                    if v[k].startswith("_"):
                        if v[k][1:] in feats: # need to strip heading _
                            v[k] = v[k][1:]
                        d_feats[k] = v[k]
            if prop_rex.match(lline):
                p = lline.split(",")
                # comment out the property
                s = "# " + p[1] + "\n\n"
            if result_rex.match(lline):
                rlist = re.split('\(|,|\)', lline)[1:-1]
                clist = list()
                for k in d_feats.keys():
                    if rlist[k] == '1':
                        clist.append(d_feats[k])
                    elif rlist[k] == '0':
                        clist.append("~"+ d_feats[k])
                s += " & ".join(clist)+"\n"
            if end_rex.match(lline):
                with open(experiment+"."+str(i), "w") as ef:
                    ef.write(s)
                    i += 1