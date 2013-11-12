#!/usr/bin/python
import sys
sys.path.append("../lib/")
import subprocess
import argparse
import random
from random import choice
from aaSeq import *
import re

sys.path.append("../lib")

def callAlign(align_version, seq1, seq2, score_file, gap):
    cmd = "%s %s %s %s %d" % (align_version, seq1, seq2, score_file, gap)

    p = subprocess.Popen(cmd.split(" "), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    
    o, e = p.communicate()
    o = o.decode()


    # returns 3 lines: alignment score then 2 sequences
    return o.rstrip().split("\n")


def homology_command(align_version, query, database, score_file, gap, p, numTrials, output_file, w):
    

    testSequences(align_version, query, database, score_file, gap, p, numTrials, output_file, w)

     
def compute_cutoff(query, database, t, p, align_version, S, gap):
    
    query = str(query)
    
    result = []
    for i in range(t):
    
        d = choice(database)
    
        l = list(d)
        random.shuffle(l)
        database_result = ''.join(l)
            
        l2 = list(query)
        random.shuffle(l2)
        query_result = ''.join(l2)    
        
        score = callAlign(align_version, query_result, database_result, S, gap)
        

        result.append(score[0])
        
    m = max(result)

    return ((1-p) * int(m))
    

def testSequences(align_version, query, database, S, gap, p, numTrials, output_file, w):

    
    query = readFA(query)[0]
    database = readFA(database)
    
    cutoff = compute_cutoff(query, database, numTrials, p, align_version, S, gap)
    print("%s: %d" % ("CUTOFF", cutoff))
    
    seed = ''
    for i in range(w):
        seed += "1"
    c_pattern = ''
    
    list_homologous = []
    for x in database:
        c_pattern = ''

        for a, b in zip(x, query):
            if str(a) == str(b):
                c_pattern += '1'
            else:
                c_pattern += '0'    

       
        if re.search(seed, c_pattern):  

            test_num = callAlign(align_version, query, x, S, gap)[0] >= cutoff
            print(test_num)
            print(callAlign(align_version, query, x, S, gap)[0])
            if test_num == True:
                list_homologous.append((x))
                    
    print("%s: %d" % ("RETAINED", len(list_homologous)))
    writeFA(list_homologous, output_file, col_width = 60)
    return list_homologous      


if __name__ == "__main__":


    align_version = "align.exe"
    parser = argparse.ArgumentParser()

    group = parser.add_mutually_exclusive_group()

    group.add_argument("-W", "--WINDOWS", action="store_true")
    group.add_argument("-X", "--OSX", action="store_true")
    group.add_argument("-L", "--LINUX", action="store_true")

    parser.add_argument('-s', '--scoring_matrix', action = 'store',
                        dest = 'score_file', default = "Blosum62.txt",
                        type = str, help = "File containing scoring matrix")

    parser.add_argument('-g', '--gap_penalty', action = 'store',
                        dest = 'gap', default = 9,
                        type = int, help = "Gap penalty (as non negative number)")    
    
    
    parser.add_argument('-o', '--output_filename', action = 'store',
                        dest = 'output_file', default = "output.fa",
                        type = str, help = "Name of output file")  
    
    parser.add_argument('-p', '--p_value', action = 'store',
                            dest = 'p', default = .05,
                            type = int, help = "P value for testing homology")
    
    parser.add_argument('-t', '--num_trials', action = 'store',
                        dest = 'num_trials', default = 30,
                        type = int, help = "Number of trials for testing homology")
    
    parser.add_argument('-w', '--length', action = 'store',
                            dest = 'w', default = 12,
                            type = int, help = "Length for required substring")    

    parser.add_argument('query', type = str, help = "Fasta file used for the query")
    parser.add_argument('database', type = str, help = "Fasta file used for the database")
    args = parser.parse_args()

    if args.WINDOWS:
        align_version = "align.exe"
    if args.OSX:
        align_version = "./align.X"
    if args.LINUX:
        align_version = "align"
    
    
    homology_command(align_version, args.query, args.database, args.score_file, args.gap, args.p, args.num_trials, args.output_file, args.w)
