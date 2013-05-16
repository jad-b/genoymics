#!/usr/bin/env python
__author__ = "Jeremy Dobbins-Bucklad"
__email__ = "jeremy.a.db@gmail.com"

import argparse
import sys
import logging
import os

parser = argparse.ArgumentParser(description=('Merges taxonomic assignments'))
# argparse opens our files through the FileType invocation
parser.add_argument('-f','--files',help='The .fna sequence file to cluster',
    nargs='+',type=argparse.FileType('r'),dest='files')
parser.add_argument('-o',help='The output directory',dest='out_dir',
    default='output')
parser.add_argument('-d','--debug',action='store_true',default=False,
    help='Print debug messages')
args = parser.parse_args()

def init(): 
    in_files = args.files   
    if not in_files:
        print "No files specified."
        return
    if len(in_files) == 1:
        print "Only one file provided; no merging possible."
        return
    
    LOG_FILENAME = 'merge_taxons.log'
    if args.debug:
        LEVEL = logging.DEBUG
    else:
        LEVEL = logging.INFO
    logging.basicConfig(filename=LOG_FILENAME,level=LEVEL,filemode='w')
    logging.info('Logging initiated')
    logging.info("Files provided: {}".format(in_files))

    return in_files


def return_best_assignment(taxon_results):
    best_assn_index = 0
    best_depth = -1
    for i,assn in enumerate(taxon_results):
        local_depth = -1
        for j,lvl in enumerate(assn):
            if len(lvl) < 4:
                break
            local_depth = j
        if local_depth > best_depth:
            best_depth = local_depth
            best_assn_index = i
    return taxon_results[best_assn_index],best_assn_index


def build_best_taxonomy(results,taxon_results):
    """
    Compare the taxonomic classifications in taxon_results
    and return the selected "best" classification.
    """
    (best_classif,index) = return_best_assignment(taxon_results)
    str = '{otu_id}\t{classif}\t{conf}\n'.format(otu_id=results[index][0],
                                                 classif=best_classif,
                                                 conf=results[index][2])
    logging.debug(str)
    return str


def main():
    in_files = init()    
    out_file=open('merged_taxons.txt','w')
    
    while True:
        results = taxon_results = []
        for i,f in enumerate(in_files):
            line = f.readline()
            # logging.debug(line)
            if line == '':
                out_file.close()
                return
            results.append(line.split())

            # Split taxonomic assignments by level
            assns = results[i][1].split(';')    
            logging.debug(assns)        
            taxon_results.append(assns)
        # logging.debug("results::{}".format(results))
        # logging.debug("Taxonomies::{}".format(taxon_results))
        #merged_line = build_best_taxonomy(results,taxon_results)
        #out_file.write(merged_line)


if __name__ == '__main__':
    main()


