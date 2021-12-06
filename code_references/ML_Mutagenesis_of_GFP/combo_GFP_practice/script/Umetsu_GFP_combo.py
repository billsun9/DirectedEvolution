import sys
import os
import numpy as np
import scipy
import combo
#import cPickle as pickle

# Note by saito:
# python Umetsu_GFP_combo.py Umetsu_GFP_BLOSUM.csv seed 2>&1 | tee Umetsu_GFP_BLOSUM.log.seed

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

# Load the data.
# X is the N x d dimensional matrix. Each row of X denotes the d-dimensional feature vector of search candidate.
# t is the N-dimensional vector that represents the corresponding negative energy of search candidates.
# ( It is of course unknown in practice. )
X, t = load_data(infile)

# Normalize the mean and standard deviation along the each column of X to 0 and 1, respectively
X = combo.misc.centering( X )

# Declare the class for calling the simulator.
# In this tutorial, we simply refer to the value of t.
# If you want to apply combo to other problems, you have to customize this class.
class simulator:
    def __init__( self, infile ):
        _, self.t = load_data(infile)
    def __call__( self, action ):
        return self.t[action]

# Declaring the policy by
policy = combo.search.discrete.policy(test_X=X)
# test_X is the set of candidates which is represented by numpy.array
# Each row vector represents the feature vector of the corresponding candidate

# set the seed parameter
print "seed: " + seed
policy.set_seed( int(seed) )

# If you want to perform the initial random search before starting the Bayesian optimization, 
# the random sampling is performed by
res = policy.random_search(max_num_probes=5, simulator=simulator(infile))
# Input
# max_num_probes: number of random search
# simulator = simulator
# output: combo.search.discreate.results (class)

# single query Bayesian search
# The single query version of COMBO is performed by 
res = policy.bayes_search(max_num_probes=95, simulator=simulator(infile), score='TS', interval=5, num_rand_basis=0) # 5000
# Input
# max_num_probes: number of searching by Bayesian optimization
# simulator: the class of simulator which is defined above
# score: the type of aquision funciton. TS, EI and PI are available
# interval: the timing for learning the hyper parameter.
#           In this case, the hyper parameter is learned at each 20 steps
#           If you set the negative value to interval, the hyper parameter learning is not performed
#           If you set zero to interval, the hyper parameter learning is performed only at the first step
# num_rand_basis: the number of basis function. If you choose 0, ordinary Gaussian process runs
#
# Note by saito:
# use num_rand_basis=0 except when computation time is a severe bottleneck.
# random feature map involves the variance of search results depending on seed parameters, 
# which makes it difficult to benckmark the performance of different feature values.


# The result of searching is summarized in the class combo.search.discrete.results.history()
# res.fx: observed negative energy at each step
# res.chosed_actions: history of choosed actions
# best_fx, best_action = res.export_all_sequence_best_fx(): 
# current best fx and current best action that has been observed until each step
# res.total_num_search: total number of search
print 'f(x)='
print res.fx[0:res.total_num_search]
best_fx, best_action = res.export_all_sequence_best_fx()
print 'current best='
print best_fx
print 'current best action='
print best_action
print 'history of chosed actions='
print res.chosed_actions[0:res.total_num_search]
# Note by saito:
# The indices of actions are 0-based.

# save the results
#res.save('test.npz')

#del res

# load the results
#res = combo.search.discrete.results.history()
#res.load('test.npz')

