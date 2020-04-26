import argparse
import numpy as np
import scipy
from scipy import stats
from numba import jit


def equivalence_class_EM():
	print("equivalence_class_EM")



@jit(nopython=True) # this is the magic numba incantation that makes things go fast
def full_model_EM(nt, label_list, ind_list, prob_list):
	num_reads = np.zeros(nt)
	eta = np.ones(nt) / float(nt) # nt is number of transcripts, this is 1/M for everything
	eta_p  = np.zeros(nt)          # this is where we will put the estimate of eta for the next iteration of the EM algorithm
	converged = False

	it = 0
	ni = len(ind_list)-1

	while not converged:
		it += 1
		for i in range(ni): # number of alignment blocks
			norm = 0.0
			for j in range(ind_list[i], ind_list[i+1]):
				P_1 = eta[label_list[j]]
				norm += P_1 * prob_list[j]

			for j in range(ind_list[i], ind_list[i+1]):
				num_reads[label_list[j]] += (eta[label_list[j]] * prob_list[j])/norm
				eta_p[label_list[j]] = num_reads[label_list[j]]/ni


		delta = np.sum(np.absolute(eta_p - eta))

		if delta < 0.001:
			converged = True
		else:
			num_reads = np.zeros(nt)

		print(it)
		eta = eta_p
		eta_p = np.zeros(nt)

	return num_reads

def calc_initial_probabilities(ref_dict, outfile):
    label_list = []
    prob_list = []
    ind_list = [0]
    mu = 200
    sd = 25
    d = scipy.stats.norm(mu, sd)

    # Read the alignment fragments and compute probabilities
    while True:
        aln_line = f.readline()
        if aln_line == "":
            break
        na = int(aln_line) # number of alignments in this alignment block
        for i in range(na):
            toks = f.readline().strip().split("\t")
            txp = toks[0]
            ori = toks[1]
            pos = int(toks[2])
            tlen, tid, eff_len = ref_dict[txp]
            label_list.append(tid)

            P_2 = float(1/eff_len)

            P_3 = d.cdf(int(tlen)-pos) if ori == "f" else d.cdf(pos + 100)

            P_4 = float(toks[3])

            prob_list.append(P_2 * P_3 * P_4)

        ind_list.append(len(prob_list))

    label_lst = np.array(label_list)
    ind_lst = np.array(ind_list)
    prob_lst = np.array(prob_list)

    print("\n\n\n\n\n=======================CALL FULL MODEL EM=========================\n\n\n\n\n")
    num_reads = full_model_EM(transcript_len, label_lst, ind_lst, prob_lst)

    keys = list(ref_dict.keys())
    effective_lens = list(ref_dict.values())
    eff_lens = [v[2] for v in effective_lens]

    with open(outfile, "w") as f1:
    	for i in range(len(num_reads)):
    		s = str(keys[i]) + "\t" + str(eff_lens[i]) + "\t" + str(num_reads[i]) + "\n"
    		f1.writelines(s)


def calc_statistics(em_file, true_counts_file):
	estimate_counts = []
	true_counts = []

	with open(em_file, "r") as f:
		for line in f:
			toks = line.strip().split("\t")
			estimate_counts.append(float(toks[2]))

	with open(true_counts_file, "r") as f:
		skip_line = f.readline()
		for line in f:
			toks = line.strip().split("\t")
			true_counts.append(float(toks[1]))

	estimate_counts = np.array(estimate_counts)
	true_counts = np.array(true_counts)

	mean_abs_error = np.mean(np.abs(_error(true_counts, estimate_counts)))
	print(mean_abs_error)

	spearman_rank_correlation_coef = scipy.stats.spearmanr(estimate_counts, true_counts)
	print(spearman_rank_correlation_coef)

	mean_absolute_percent_error = np.mean(np.abs(_percentage_error(true_counts, estimate_counts)))
	print(mean_absolute_percent_error)

def _error(actual: np.ndarray, predicted: np.ndarray):
    """ Simple error """
    return actual - predicted

def _percentage_error(actual: np.ndarray, predicted: np.ndarray):
    """
    Percentage error

    Note: result is NOT multiplied by 100
    """
    EPSILON = 1e-10
    return _error(actual, predicted) / (actual + EPSILON)

def memoizer(f):
        memolist = [-1] * 1000
        def _f(l):
                if memolist[l] == -1:
                        memolist[l] = f(l)
                return memolist[l]
        return _f

@memoizer
def get_eff_len_short(l):
        mu = 200
        sd = 25
        d = scipy.stats.norm(mu, sd)
        # the (discrete) distribution up to 200
        p = np.array([d.pdf(i) for i in range(l+1)])
        # re-normalize so this is a proper distribution
        p /= p.sum()
        # the expected value of a distribution f is \sum_{i=0}^{max_i} f(i) * i
        cond_mean = np.sum([i * p[i] for i in range(len(p))])
        return l - cond_mean

# the "main" function, this is the one you call to get the effective length
# efficiently
def get_eff_len(l):
        mu = 200
        return l - mu if l >= 1000 else get_eff_len_short(l)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='code')
    parser.add_argument('--i', type=str, default='data/alignments.cs423',
                        help='path to format alignment file')

    parser.add_argument('--out', type=str, default='quant.tsv',
                        help='path to file where output should be written')

    parser.add_argument('--eqc', type=str, default=False,
                        help='wether to perform an equivalence class EM or a full-model EM')

    parser.add_argument('--stats', type=str, default=False,
                        help='True if you want to perform statistics on EM predictions vs actual data')

    parser.add_argument('--em_estimation_file', type=str, default='quant.tsv',
                        help='The output tsv file genereated using the EM Model to predict counts and effective length')

    parser.add_argument('--true_counts_file', type=str, default='data/true_counts.tsv',
                        help='The file that holds the true counts of each transcript')


    args = parser.parse_args()

    if args.stats:
    	calc_statistics(args.em_estimation_file, args.true_counts_file)
    else: 
	    ref_dict = {}

	    # Read the transcripts and save in ref_dict
	    f = open(args.i, "r")
	    transcript_len = int(f.readline())

	    for i in range(transcript_len):
	    	toks = f.readline().strip().split("\t")
	    	# transcript name => (actual length, index, effective length)
	    	ref_dict[toks[0]] = [toks[1], i, get_eff_len(int(toks[1]))]

	    print("CALC INITIAL PROBS")
	    calc_initial_probabilities(ref_dict, args.out)

	    # equivalence_class_EM() if args.eqc else full_model_EM()

