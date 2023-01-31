'''
    Romain Derelle (Imperial College London)
    contact: romain.derelle@gmail.com
'''

import sys
import os
import csv
import collections
import math
import statistics
import kneed
import argparse



def parse_args():
    # define and parse command-line arguments
    parser = argparse.ArgumentParser(description='       SCBD-LCBD analysis', add_help=False, formatter_class=argparse.RawTextHelpFormatter, epilog=' \n')
    common = parser.add_argument_group()
    common.add_argument('-m',         help='path to the normalised matrix file [required]', metavar='')    
    common.add_argument('-h','-help', action="help", help="show this help message and exit")
    args = parser.parse_args()
    
    return args.m


file_matrix = parse_args()

## sanity check
if not file_matrix:
    sys.exit('\n            ERROR: please specify a matrix file (-m option)\n\n')
elif not os.path.isfile(file_matrix):
    sys.exit('\n            ERROR: the path to the matrix file is incorrect\n\n')


## get all data from CSV file
print('read input matrix')
nb_line = 0
counter_feature  = 0
feature_names  = dict()
feature_values = dict()
raw_data = dict()

file_content = csv.reader(open(file_matrix), delimiter=',')
for line in file_content: 
    nb_line += 1
    # extract column names
    if nb_line == 1:
        sample_names = dict()
        total = dict()
        for i,k in enumerate(line):
            sample_names[i] = k 
            total[sample_names[i]] = 0
    
    # extract feature values
    elif nb_line > 1:
        counter_feature += 1
        name = line[0]
        feature_names[counter_feature] = name
        # get data
        raw_data[counter_feature] = line[1:]
        feature_values[counter_feature]   = dict()
        for n in range(1, len(line)):
            nb_reads = line[n]
            sample = sample_names[n]
            feature_values[counter_feature][sample] = float(nb_reads)
            total[sample] += float(nb_reads)
            


## perform Hellinger transformation (square root of the fraction)
print('perform Hellinger transformation')
hellinger = collections.defaultdict(dict)
for feature in feature_values:
    for sample in feature_values[feature]:
        nn = math.sqrt(feature_values[feature][sample] / total[sample])
        hellinger[feature][sample] = nn
        
        
## get mean value for each feature
print('calculate means and SS value')
feature_means = dict()
for feature in hellinger:
    feature_means[feature] = statistics.mean([v for v in hellinger[feature].values()])


## build S matrix
S_matrix = collections.defaultdict(dict) 
for feature in hellinger:
    for sample in hellinger[feature]:
        S_matrix[sample][feature] = math.pow(hellinger[feature][sample] - feature_means[feature], 2)    ## invert rows and columns


## calculate SS value
SS_total = 0
for sample in S_matrix:
    for feature in S_matrix[sample]:
       SS_total += S_matrix[sample][feature]


## calculate SCBD value for each feature
print('calculate SCBD values')
SCBD = dict()
for feature in feature_names:
    tmp_sum = 0
    for sample in S_matrix:
        tmp_sum += S_matrix[sample][feature]
    SCBD[feature] = tmp_sum / SS_total


## sort the SCBD by values (decreasing order)
l_sorted_SCBD = list()
for key, value in sorted(SCBD.items(), key=lambda x: x[1], reverse=True):
    l_sorted_SCBD.append(key)


## save SCBD values
print('-> save SCBD values to file')
out_SCBD = open('out_SCBD_values.tsv', 'w+')
for feature in l_sorted_SCBD:
    out_SCBD.write(feature_names[feature] + '	' + str(SCBD[feature]) + '\n')
out_SCBD.close()


## save sorted matrix
print('-> save sorted matrix to file')
out_matrix = open('out_SCBD_matrix.csv', 'w+')
for feature in l_sorted_SCBD:
    out_matrix.write(feature_names[feature] + ',' + ','.join(raw_data[feature]) + '\n')
out_matrix.close()    





