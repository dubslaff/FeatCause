import pandas as pd

experiments = [
    [['CFDP', 'L'], 'cfdp', 'cfdpNegated',],
    [['Elevator\_1', 'L'], 'elevator', 'elevatorNegated',],
    [['Minepump', 'L'], 'minepump', 'minepumpNegated',],
    [['BSN', 'L'], 'BSN', 'BSNNegated',],
    [['Aircraft', 'L'], 'merged_single', 'merged_singleNegated',],
    [['lrzip', 'T'], 'lrzipPaper', 'lrzipPaperNegated'],
    [['x264', 'T'], 'x264Paper', 'x264PaperNegated'],
    [['DUNE', 'T'], 'dunePaper', 'dunePaperNegated'],
    [['BerkeleyDB','T'], 'berkeleyDBCPaper', 'berkeleyDBCPaperNegated'],
    [['LLVM', 'T'], 'llvmPaper', 'llvmPaperNegated'],
    [['Apache', 'P'], 'Apache_PFW', 'Apache_PPW', 'Apache_PFWNegated', 'Apache_PPWNegated',],
    [['Curl', 'M'], 'Curl_MEMFW', 'Curl_MEMPW', 'Curl_MEMFWNegated', 'Curl_MEMPWNegated',],
    [['h264','PM'], 'h264_MEMFW', 'h264_MEMPW', 'h264_PFW', 'h264_PHO', 'h264_PHS', 'h264_PPW',
      'h264_MEMFWNegated', 'h264_MEMPWNegated', 'h264_PFWNegated', 'h264_PHONegated', 'h264_PHSNegated', 'h264_PPWNegated',],
    [['LinkedList', 'B'], 'LinkedList_BSDW', 'LinkedList_BSFW', 'LinkedList_BSPW',
      'LinkedList_BSDWNegated', 'LinkedList_BSFWNegated', 'LinkedList_BSPWNegated',],
    [['Linux', 'B'], 'Linux_BSDW', 'Linux_BSFW', 'Linux_BSPW',],
    [['PKJab', 'B'], 'PKJab_BSFW', 'PKJab_BSPW', 'PKJab_BSFWNegated', 'PKJab_BSPWNegated',],
    [['Prevaylar', 'B'], 'Prevaylar_BSFW', 'Prevaylar_BSDW', 'Prevaylar_BSPW',
      'Prevaylar_BSFWNegated', 'Prevaylar_BSDWNegated', 'Prevaylar_BSPWNegated',],
    [['SNW', 'B'], 'SNW_BSFW', 'SNW_BSDW', 'SNW_BSPW', 'SNW_BSFWNegated', 'SNW_BSDWNegated', 'SNW_BSPWNegated',],
    [['SQLite', 'M'], 'sqlite_MEMFW', 'sqlite_MEMPW', 'sqlite_MEMPWNegated',],
    [['WGet', 'M'], 'wget_MEMFW', 'wget_MEMPW', 'wget_MEMFWNegated', 'wget_MEMPWNegated',],
    [['ZipMe', 'BP'], 'ZipMe_BSFW', 'ZipMe_BSDW', 'ZipMe_BSPW', 'ZipMe_PFW', 'ZipMe_PFami', 'ZipMe_PPW',
      'ZipMe_BSFWNegated', 'ZipMe_BSDWNegated', 'ZipMe_BSPWNegated', 'ZipMe_PFWNegated', 'ZipMe_PFamiNegated', 'ZipMe_PPWNegated',],
    [['Elevator\_2', 'P'], 'Elevator_PFW', 'Elevator_PPW', 'Elevator_PFWNegated', 'Elevator_PPWNegated',],
    [['E-Mail', 'P'], 'EMail_PFW', 'EMail_PPW', 'EMail_PFWNegated', 'EMail_PPWNegated',],
]

experiment_dir = '../tmp/statistic_files/'


def construct_nfp_entry(nfps):
    nfp = ''
    if 'P' in nfps:
        nfp += 'P'
    if 'B' in nfps:
        nfp += 'B'
    if 'M' in nfps:
        nfp += 'M'
    if 'T' in nfps:
        nfp += '$\\vartheta_P$'
    if 'L' in nfps:
        nfp += '$\\varphi$'
    return nfp


def create_average():
    df = {'system': [],
          '# experiments': [],
          'V': [],
          'F': [],
          'NFP': [],
          'time': [],
          'avg |E|': [],
          'avg At': [],
          'avg MaxCov': [],
          'N->K': [],
         }
    for exp in experiments:
        frames = []
        for e in exp[1:]:
            dataframe_exp = pd.read_csv(experiment_dir + e + 'Stats.csv')
            if 'n' not in dataframe_exp or 'k' not in dataframe_exp:
                print(exp)
            frames.append(dataframe_exp)
        dataframe_exp = pd.concat(frames)
        if 'TimeAtomicMinimals' not in dataframe_exp:
            print(exp, e)
        df['system'].append(exp[0][0])
        df['# experiments'].append(len(dataframe_exp.index))
        df['V'].append(dataframe_exp['NumberValidConfigurations'].mean())
        df['F'].append(dataframe_exp['NumberFeatures'].mean())
        df['NFP'].append(construct_nfp_entry(exp[0][1]))
        df['time'].append(dataframe_exp['TimeAtomicMinimals'].sum())
        df['avg |E|'].append(dataframe_exp['EffectSetSize'].mean())
        df['avg At'].append(dataframe_exp['NumberAtomics'].mean())
        df['avg MaxCov'].append(dataframe_exp['NumberMinimals'].mean())
        df['N->K'].append(1 - (dataframe_exp['k'].mean()/dataframe_exp['n'].mean()))
    df = pd.DataFrame.from_dict(df)
    print(df)
    return df


def main():
    df = create_average()
    df['V'] = df['V'].astype(int)
    df['F'] = df['F'].astype(int)
    df.to_latex('statisticsTable.tex', index=False, header=['','\#', '$|\\ValidFeat|$', '$|\\Feat|$', 'Property', 'Time', '$\\Effects$', '$\\Causes$', '$\\mCauses$', 'DLS' ], float_format='{%.2f}', column_format='l|rrrlr|rrrr', escape=False)
    pass


def parse_arguments():
    import argparse
    ap = argparse.ArgumentParser(description="A script to generate the table in the evaluation of the ICSE'22 paper.
                                              Note: First, run the evaluation_runner script, which creates all statistics regarding the experiments,
                                              then you can build the evaluation table.")


if __name__ == '__main__':
    args = parse_arguments()
    main()
