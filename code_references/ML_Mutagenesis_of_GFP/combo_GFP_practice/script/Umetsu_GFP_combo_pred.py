import sys
import os
import numpy as np
import scipy
import combo
#import cPickle as pickle

# Note by saito:
# python Umetsu_GFP_combo_pred.py Umetsu_GFP_BLOSUM_pred.csv seed npred 2>&1 | tee Umetsu_GFP_BLOSUM_pred.log.seed.npred

ntrain = 155

# set the stdout buffer size to 0, which is useful for "| less" or "| tee log" 
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)


def load_data(infile):
    print "loading " + infile
    with open(infile) as f:
        ncol = len(f.readline().split(','))
    print "ncol " + str(ncol)
    A = np.asarray( np.loadtxt(infile,skiprows=1,usecols=range(1,ncol),delimiter=',') )
    X = A[:,0:ncol-2]
    t = A[:,ncol-2]
    return X, t
# Note by saito:
# The values of t must be negative. 
# The algorithm searches for the action which has the largest value of t 
# i.e. the minimum absolute value of t

argvs = sys.argv
infile = argvs[1]
seed = argvs[2]

npred = int(argvs[3])

# Load the data.
# X is the N x d dimensional matrix. Each row of X denotes the d-dimensional feature vector of search candidate.
# t is the N-dimensional vector that represents the corresponding negative energy of search candidates.
# ( It is of course unknown in practice. )
X, t = load_data(infile)

# Normalize the mean and standard deviation along the each column of X to 0 and 1, respectively
X = combo.misc.centering( X )

train_X = X[0:ntrain, :]
train_t = t[0:ntrain]
test_X = X[ntrain:, :]
test_t = t[ntrain:] # dummy

# Declaring the policy by
policy = combo.search.discrete.policy(test_X=test_X)
# test_X is the set of candidates which is represented by numpy.array
# Each row vector represents the feature vector of the corresponding candidate

# set the seed parameter
print "seed: " + seed
policy.set_seed( int(seed) )

training = combo.variable(X=train_X, t=train_t)

actions, ff = policy.bayes_search(training=training, max_num_probes=1, num_search_each_probe=npred, simulator=None, score='TS', interval=0, num_rand_basis = 5000) # 5000
# Note by saito:
# In this problem, computation time is a severe bottleneck due to the large size of test data (~20^4). 
# Therefore, we use num_rand_basis=5000 instead of num_rand_basis=5000. 

if npred == 1:
    print "ff:"
    for i in xrange(len(ff)):
        print str(i) + "\t" + str(ff[i])
    print "end"

print "action:"
print actions
print "end"
