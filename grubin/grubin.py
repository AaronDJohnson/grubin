import numpy as np
import pandas as pd
import os
from glob import glob


def load_chain(chain_loc, burn=0, thin=1, numpy=False):
    """
    Import chains to pandas.

    Input:
        chaindir (string): folder path of chain to be loaded
        chain_loc (string): file path of chain to be loaded
        burn (float) [0]: percent of chain to drop from the beginning
        thin (float) [1]: take every (#)rd sample -- 2 = every other sample
        numpy (bool) [False]: return as numpy array?

    Output:
        chain_df (DataFrame): chain output as pandas dataframe
    """
    # import params for headers
    chaindir = os.path.dirname(chain_loc)
    param_file = chaindir + '/pars.txt'
    with open(param_file, 'r') as f:
        params = [line.rstrip('\n') for line in f]
    params += ['lnpost', 'lnlike', 'accept_rate', 'pt_accept']

    # import data
    chain_raw = pd.read_csv(chain_loc, sep='\t', dtype=float, header=None).values
    chain_df = pd.DataFrame(data=chain_raw, columns=params)
    if numpy:
        chain_df = chain_df.to_numpy()
    return chain_df


def grubin(data, M=4, burn=0.25, threshold=1.1):
    """
    Gelman-Rubin split R hat statistic to verify convergence.

    See section 3.1 of https://arxiv.org/pdf/1903.08008.pdf.
    Values > 1.1 => recommend continuing sampling due to poor convergence.

    Input:
        data (ndarray): consists of entire chain file
        pars (list): list of parameters for each column
        M (integer): number of times to split the chain
        burn (float): percent of chain to cut for burn-in
        threshold (float): Rhat value to tell when chains are good

    Output:
        Rhat (ndarray): array of values for each index
        idx (ndarray): array of indices that are not sampled enough (Rhat > threshold)
    """
    burn = int(burn * data.shape[0])  # cut off burn-in

    try:
        data_split = np.split(data[burn:,:-2], M)  # cut off last two columns
    except:
        # this section is to make everything divide evenly into M arrays
        P = int(np.floor((len(data[:, 0]) - burn) / M))  # nearest integer to division
        X = len(data[:, 0]) - burn - M * P # number of additional burn in points
        burn += X  # burn in to the nearest divisor
        burn = int(burn)

        data_split = np.split(data[burn:,:-2], M)  # cut off last two columns

    N = len(data[burn:, 0])
    data = np.array(data_split)

    # print(data_split.shape)

    theta_bar_dotm = np.mean(data, axis=1)  # mean of each subchain
    theta_bar_dotdot = np.mean(theta_bar_dotm, axis=0)  # mean of between chains
    B = N / (M - 1) * np.sum((theta_bar_dotm - theta_bar_dotdot)**2, axis=0)  # between chains

    # do some clever broadcasting:
    sm_sq = 1 / (N - 1) * np.sum((data - theta_bar_dotm[:, None, :])**2, axis=1)
    W = 1 / M * np.sum(sm_sq, axis=0)  # within chains
    
    var_post = (N - 1) / N * W + 1 / N * B
    Rhat = np.sqrt(var_post / W)

    idx = np.where(Rhat > threshold)[0]  # where Rhat > threshold

    return Rhat, idx


def check_chains(chaindir=None):
    """
    Searches through subdirectories for chain files.
    """
    if not chaindir:
        fdir = os.getcwd()
    else:
        fdir = chaindir
    folder_names = [x[0] for x in os.walk(fdir)]
    if not folder_names:
        print('No folders found in given directory. Does the path exist?')
    for folder in folder_names:
        files = sorted(glob(folder + '/chain*'))
        if not files:
            print('No chain files found in', folder)
        else:
            done_flag = True
            for fname in files:
                with open(fname, 'r') as f:
                    try:
                        data = load_chain(fname, numpy=True)
                    except:
                        print('File', fname, 'is empty. Continuing...')
                        continue
                rhat, idx = grubin(data)
                if idx.size > 0:
                    done_flag = False
                    print('File', fname, 'has not converged yet.')
                else:
                    print('File', fname, 'has converged.')
            if done_flag:
                print('All files under the Gelman-Rubin statistic threshold.')
            else:
                print('See above files that need more sampling.')

    

if __name__ == '__main__':
    find_files()









