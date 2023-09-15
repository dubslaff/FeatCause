experiments = [
    'RWSPL/llvm/llvmPaper',
    'RWSPL/lrzip/lrzipPaper',
    'RWSPL/x264/x264Paper',
    'RWSPL/Dune/dunePaper',
    'RWSPL/berkeleyDBC/berkeleyDBCPaper',
    'PCSimp/E1/Apache_P/Apache_PFW',
    'PCSimp/E1/Apache_P/Apache_PPW',
    'PCSimp/E1/Curl_MEM/Curl_MEMFW',
    'PCSimp/E1/Curl_MEM/Curl_MEMPW',
    'PCSimp/E1/Elevator_P/Elevator_PFW',
    'PCSimp/E1/Elevator_P/Elevator_PPW',
    'PCSimp/E1/EMail_P/EMail_PFW',
    'PCSimp/E1/EMail_P/EMail_PPW',
    'PCSimp/E1/h264_MEM/h264_MEMFW',
    'PCSimp/E1/h264_MEM/h264_MEMPW',
    'PCSimp/E1/h264_P/h264_PFW',
    'PCSimp/E1/h264_P/h264_PHO',
    'PCSimp/E1/h264_P/h264_PHS',
    'PCSimp/E1/h264_P/h264_PPW',
    'PCSimp/E1/LinkedList_BS/LinkedList_BSDW',
    'PCSimp/E1/LinkedList_BS/LinkedList_BSFW',
    'PCSimp/E1/LinkedList_BS/LinkedList_BSPW',
    'PCSimp/E1/Linux_BS/Linux_BSDW',
    'PCSimp/E1/Linux_BS/Linux_BSFW',
    'PCSimp/E1/Linux_BS/Linux_BSPW',
    'PCSimp/E1/llvm_MEM/llvm_MEMFW',
    'PCSimp/E1/llvm_MEM/llvm_MEMPW',
    'PCSimp/E1/llvm_P/llvm_PFW',
    'PCSimp/E1/llvm_P/llvm_PHO',
    'PCSimp/E1/llvm_P/llvm_PHS',
    'PCSimp/E1/llvm_P/llvm_PPW',
    'PCSimp/E1/PKJab_BS/PKJab_BSFW',
    'PCSimp/E1/PKJab_BS/PKJab_BSPW',
    'PCSimp/E1/Prevaylar_BS/Prevaylar_BSFW',
    'PCSimp/E1/Prevaylar_BS/Prevaylar_BSDW',
    'PCSimp/E1/Prevaylar_BS/Prevaylar_BSPW',
    'PCSimp/E1/SNW_BS/SNW_BSFW',
    'PCSimp/E1/SNW_BS/SNW_BSDW',
    'PCSimp/E1/SNW_BS/SNW_BSPW',
    'PCSimp/E1/sqlite_MEM/sqlite_MEMFW',
    'PCSimp/E1/sqlite_MEM/sqlite_MEMPW',
    'PCSimp/E1/wget_MEM/wget_MEMFW',
    'PCSimp/E1/wget_MEM/wget_MEMPW',
    'PCSimp/E1/ZipMe_BS/ZipMe_BSFW',
    'PCSimp/E1/ZipMe_BS/ZipMe_BSDW',
    'PCSimp/E1/ZipMe_BS/ZipMe_BSPW',
    'PCSimp/E1/ZipMe_P/ZipMe_PFW',
    'PCSimp/E1/ZipMe_P/ZipMe_PFami',
    'PCSimp/E1/ZipMe_P/ZipMe_PPW',
    'ProVeLines/cfdp/cfdp',
    'ProVeLines/elevator/elevator',
    'ProVeLines/minepump/minepump',
    'Prism/BSN/BSN',
    'Prism/aircraft/merged_single',
]

experiment_dir = '../benchmark/examples/'


class bc:
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

    def color_string(color, string):
        return color + str(string) + bc.ENDC


def check_experiment_list():
    """
    This method checks if all listed experiments do actually exist (at the given path)
    In the case an experiment is not found, the experiment (with the given path) is printed.
    At the end, a counter is printed how many experiments do not exist.
    """
    import os
    non_existing_exp = 0
    for file_name in experiments:
        if not os.path.isfile(experiment_dir + file_name + '.0'):
            print(bc.color_string(bc.FAIL, experiment_dir + file_name), 'is not an existing experiment.')
            non_existing_exp += 1
    if non_existing_exp == 0:
        print(bc.color_string(bc.OKGREEN, 'All experiments exist.'))
    elif non_existing_exp == 1:
        print(bc.color_string(bc.WARNING, non_existing_exp), 'experiment does not exist.')
    else:
        print(bc.color_string(bc.WARNING, non_existing_exp), 'experiments do not exist.')


def worker():
	import subprocess
	while True:
		item = q.get()
		print('Starting', bc.color_string(bc.BOLD, item[1]))
		subprocess.call(item[0], shell=True)
		q.task_done()
		print(bc.color_string(bc.OKGREEN, item[1]), 'done.')


def worker_fscore():
    import subprocess
    while True:
        item = q_fscore.get()
        print('Starting', bc.color_string(bc.BOLD, item[1]))
        subprocess.call(item, shell=True)
        q_fscore.task_done()
        print(bc.color_string(bc.OKGREEN, item[1]), 'done.')


def run_all_experiments(multi_threaded):
    import queue
    import threading
    import multiprocessing

    global q
    q = queue.Queue()
    for file_name in experiments:
        experiment = experiment_dir + file_name
        run_string = ['python3 main.py ' + experiment + ' -stat -m -d > /dev/null', experiment]
        q.put(run_string)
        run_string = ['python3 main.py ' + experiment + ' -stat -m -d -n > /dev/null', experiment]
        q.put(run_string)
    if multi_threaded:
        cpus = multiprocessing.cpu_count()
    else:
        cpus = 1
    print("Creating %d threads." % cpus)
    for i in range(cpus):
        t = threading.Thread(target=worker)
        t.daemon = True
        t.start()
    q.join()


def run_timing_experiments(multi_threaded):
    import queue
    import threading
    import multiprocessing

    global q
    q = queue.Queue()
    for file_name in experiments:
        experiment = experiment_dir + file_name
        run_string = ['python3 main.py ' + experiment + ' -tstat -na -v > /dev/null', experiment]
        q.put(run_string)
        if 'Linux_BSDW' not in experiment and 'sqlite_MEMFW' not in experiment:
            run_string = ['python3 main.py ' + experiment + ' -tstat -na -v -n > /dev/null', experiment]
            q.put(run_string)
    if multi_threaded:
        cpus = multiprocessing.cpu_count()
    else:
        cpus = 1
    for i in range(cpus):
        t = threading.Thread(target=worker)
        t.daemon = True
        t.start()
    q.join()


def run_fscore_experiments(multi_treaded):
    import queue
    import threading
    import multiprocessing
    import subprocess

    # get minepump experiment
    experiment = experiment_dir + 'ProVeLines/minepump/minepump'
    
    inc_number = 500
    global q_fscore
    q_fscore = queue.Queue()
    
    if multi_treaded:
        cpus = multiprocessing.cpu_count()
        for i in range(41):
            run_string = ['python3 main.py -m -s ' + str(i) + ' -inc ' + str(inc_number) + ' ' + experiment + '  > /dev/null', experiment]
            q_fscore.put(run_string)
        for i in range(cpus):
            t = threading.Thread(target=worker_fscore)
            t.daemon = True
            t.start()
        q_fscore.join()
    else:
        run_string = 'python3 main.py -m -inc ' + str(inc_number) + ' ' + experiment + '  > /dev/null'
        subprocess.call(run_string, shell=True)


def main():
    args = parse_arguments()
    if args.check:
        check_experiment_list()
    if args.run:
        run_all_experiments(args.multi_threading)
        run_fscore_experiments(args.multi_threading)
    if args.time_run:
        run_timing_experiments(args.multi_threading)
    if args.fscore and not args.run:
        run_fscore_experiments(args.multi_threading)


def parse_arguments():
    import argparse
    ap = argparse.ArgumentParser(description='A script to execute all experiments')
    ap.add_argument('-c', '--check', action='store_true', help='Checks if all listed experiments exist.')
    ap.add_argument('-r', '--run', action='store_true', help='Use this flag to run all experiments.')
    ap.add_argument('-t', '--time_run', action='store_true', help='Use this flag to run all experiments and get the run time statistics. Note: You should not use multi_threading to reduce noise on the timing data.')
    ap.add_argument('-m', '--multi_threading', action='store_true', help='Run all experiments multithreaded.')
    ap.add_argument('-f', '--fscore', action='store_true', help='Run the fscore experiment.')
    args = ap.parse_args()
    return args


if __name__ == '__main__':
	main()
