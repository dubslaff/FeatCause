# a simple hack to run ProVeLines and generate DNF effects that FeatCause understands
# for strange reasons it does only work without explicit providing the provelines dictionary

import argparse
import re
import shlex
import subprocess

property_rex = re.compile('Checking LTL property.*')


# run ProVeLines
def run_provelines(fmodel, model, ltl, logf):
    with open(logf, "w") as log:
        command_line = "../src/provelines"+" "+ \
        "-check"+" "+ "-nt"+" "+ "-exhaustive"+" "+ "-fm"+" "+ str(fmodel)+" "+ "-ltl"+" "+ "'" + str(ltl) + "'"+" "+ str(model)
        arguments = shlex.split(command_line)
        subprocess.call(arguments, stdout=log, stderr=None)


# interpret result of ProVeLines
def extract_log(logf, outf):
    with open(logf, "r") as log:
        with open(outf, "w") as out:
            for lline in log:
                if not lline.startswith("   ") and not ('Checking' in lline):
                    continue
                if property_rex.match(lline):
                    lline = lline.replace("..", "")
                    # comment out the property
                    lline = "# " + lline.replace("Checking LTL property ", "") + "\n"
                else:
                    lline = lline.replace("All the products", "1")
                    lline = lline.replace("   ", "")
                    lline = lline.replace("!", "~")  # different format of negation
                out.write(lline)


# main program
if __name__ == '__main__':
    # define input arguments for the extractor script
    argparser = argparse.ArgumentParser(description="PVL Extract - generate DNF effects for FeatCause")
    argparser.add_argument("experiment", help="path to feature model (.tvl) and behavioral model (.pml), e.g., ./advanced/elevator")
    argparser.add_argument("props", help="list of LTL properties (.ltl)")
    args = argparser.parse_args()
    filename_experiment = args.experiment
    filename_props = args.props

    # read properties and iterate through them
    with open(filename_props) as f:
        i = 0
        for line in f:
            filename_log = filename_experiment + "." + str(i) + ".log"
            line = line.replace("\n", "")
            if line.startswith("#") or len(line) == 0 or line == "" or line == " ":
                continue
            print(str(i) + ": checking " + line)
            run_provelines(filename_experiment+".tvl", filename_experiment+".pml", line, filename_log)
            extract_log(filename_log, filename_experiment+"."+str(i))
            i += 1