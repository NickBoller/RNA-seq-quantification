import argparse
import numpy as np
import scipy
from scipy import stats
from numba import jit
@jit(nopython=True) # this is the magic numba incantation that makes things go fast


def equivalence_class_EM():
	print("equivalence_class_EM")




def full_model_EM(nt, label_list, ind_list, probs_list):
	eta = np.ones(nt) / float(nt) # nt is number of transcripts, this is 1/M for everything
	eta_p  = np.zeros(nt)          # this is where we will put the estimate of eta for the next iteration of the EM algorithm
	converged = False
	print(eta)
	# Calculate P1

	it = 0
	ni = len(ind_list)-1

	while True:
		it += 1
		for i in range(ni): # number of alignment blocks
			norm = 0.0
			for j in range(ind_list[i], ind_list[i+1]):
				# here, j is the index we can use to index *directly* into label_list and prob_list
				# iterate over this one time to compute the normalizer / denominator for each alignment block
				for j in range(ind_list[i], ind_list[i+1]):
				# here, j is the index we can use to index *directly* into label_list and prob_list
				# iterate over this a second time to proportionally allocate this fragment's contribution to each transcript
				# to which it aligns
		# Now, outside of this loop, we've done an iteration of the E step, so we can compute
		# our new eta and check for convergence etc.  (if converged, then break!)
		# at the end of this loop, before returning to the top, set eta = eta_p and clear out eta_p to 0


def calc_initial_probabilities(ref_dict):
	# assume we've already read in the "header" of the 
    # alignment file containing the number of transcripts and the 
    # name and length of each.  At this point, ref_dict  is a 
    # dictionary from transcript name to a tuple of (transcript length, transcript index) 

    '''
		- P1: ith entry of the eta array (at first it is set to 1/M)
		- P2: 1/(effective length of transcript) -> 1/get_eff_length()
		- P3: 
		- P4: given -> 4th column of alignment record

    '''

    label_list = []
    prob_list = []
    ind_list = [0]

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
            tlen, tid = ref_dict[txp]
            label_list.append(tid)

            P_2 = float(1/get_eff_len(int(tlen)))
            P_3 = scipy.stats.norm.cdf(int(tlen)-pos)
            P_4 = float(toks[3])

            prob_list.append(P_2 * P_3 * P_4)

        ind_list.append(len(prob_list))

        label_list = np.array(label_list)
        ind_list = np.array(ind_list)
        prob_list = np.array(prob_list)
        print("label list: ")
        print(label_list)
        print("probs list: ")
        print(prob_list)
        print("ind list: ")
        print(ind_list)

        full_model_EM(transcript_len, label_list, ind_list, prob_list)



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
    parser.add_argument('--inn', type=str, required=True,
                        help='path to format alignment file')

    parser.add_argument('--out', type=str, required=True,
                        help='path to file where output should be written')

    parser.add_argument('--eqc', type=str, default=False,
                        help='wether to perform an equivalence class EM or a full-model EM')

    args = parser.parse_args()


    ref_dict = {}

    # Read the transcripts and save in ref_dict
    f = open(args.inn, "r")
    transcript_len = int(f.readline())
    print("M = " + str(transcript_len))

    for i in range(transcript_len):
    	toks = f.readline().strip().split("\t")
    	ref_dict[toks[0]] = (toks[1], i)

    print("ref dict: " + str(ref_dict))

    calc_initial_probabilities(ref_dict)

    # equivalence_class_EM() if args.eqc else full_model_EM()

