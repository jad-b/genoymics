#!/usr/bin/env python
__author__ = "Jeremy Dobbins-Bucklad"
__email__ = "jeremy.a.db@gmail.com"

"""
master_hiClust.py coordinates the actions, and is responsible for calling each step. 
This is to provide easy high-level abstraction, as well as to hopefully make each 
piece modular.
"""

# Standard Python modules
import argparse 
import subprocess           
import re
import glob
import os
import logging


#Setup cmd line parser
parser = argparse.ArgumentParser(
    description=('Runs hierarchical clustering on a given sequence file'))
parser.add_argument('sequences',
                    help='The .fna sequence file to cluster')
parser.add_argument('processors',type=int,
                    help='The number of processes to run via MPI')
parser.add_argument('-o','--output_directory',help='The top-level output directory',
                    default='output',dest='output_directory')
parser.add_argument('-1','--first',default='0.75',type=float,
                    help='Similarity to cluster for first pass')
parser.add_argument('-2','--second',default='0.97',type=float,
                    help='Similarity to cluster for second pass')
parser.add_argument('-d','--debug',action='store_true',default=False,
                    help='Print debug messages')
parser.add_argument('--skipUclust',action='store_true',default=False,
                    help='Skips first clustering. Only useful for testing')
parser.add_argument('--keep',action='store_true',default=False,dest='keep',
                    help='Deletes intermediate files from this run')
parser.add_argument('-al','-append_log',action='store_true',default=False,
                    dest='store_log')
args = parser.parse_args()

if args.debug:
    LEVEL = logging.DEBUG
else:
    LEVEL = logging.INFO

if args.store_log:
    FILEMODE='a'
    FILENAME='hc.log' 
else:
    FILEMODE='w'
    #log_files=glob.glob('hc_*.log')
    #FILENAME='hc_{}.log'.format(len(log_files))
    FILENAME='hc.log'

logging.basicConfig(filename=FILENAME,
                    filemode=FILEMODE,
                    format='%(asctime)s : %(levelname)-8s: %(name)s: %(message)s',
                    level=LEVEL,
                    datefmt='%Y %b %d, %H:%M:%S')
logger = logging.getLogger('HiClust')
logger.info('logger initiated')

# variables
sim_first = args.first
first_num = str(int(sim_first*100))
sim_second = args.second
second_num = str(int(sim_second*100))
num_seqs = 0    # Number of sequences provided
method = 'pick_otus.py'
debug = args.debug


class Divvy:
    def __init__(self):
        self.logger = logging.getLogger('Divvy')
    
    def validate_directory(self,dir):
        """
        Check if 'dir' is a valid directory.
        If it exists, return it.
        If it doesn't exist, attempt to make it.
        If an error occurs, raise an exception.
        """
        try:        
            if os.path.exists(dir):
                return dir
            else:
                os.makedirs(dir)
                return dir
        except OSError as e:
            logger.error('Problem with creation of {}'.format(dir))
            raise 
    
    def retrieve_seq_by_id(self,seq_id,seq_file):
        # Iterate through, looking for seq_id
        for line in seq_file:
            if line[0] != '>': # Skip non-headers for efficiency
                continue
            if seq_id in line:
                break
    
        if line == '': # Check for EOF
            return None
    
        target_seq = line           
        for line in seq_file:
            if line[0] == '>': # Quit on a header
                break       
            target_seq += line 
    
        seq_file.seek(0) # Reset file position
        return target_seq
    
    def divvy_otus(self,divvy_dir, otu_map, seq_file, seqs_per_file ):
        """
        - divvy_dir is the output directory for divvied sequences
        - otu_dir is the directory containing the otu_map
        - seq_file is the base file of sequences
        - seqs_per_file is the target number of sequences for each file
        Reads the assigned sequence IDs for each OTU, and copies them from
        the sequence file to a new file. Attempts to only read a certain number
        of sequences per file, but will not split an OTU.
        """
        seqs_dir,seqs_filename = os.path.split(seq_file)
        seqs_name,seq_ext = os.path.splitext(seqs_filename)
        base_otu_id = None
        curr_otu_id = 0
        seqs_read = 0
        file_count = 0
        divvied_seqs = ""
    
        with open(otu_map,'r') as otus,open(seq_file,'r') as seqs:
            for otu_line in otus:
                otu_info = otu_line.rstrip('\n').split('\t')
                curr_otu_id = otu_info[0]
                if base_otu_id is None: 
                    base_otu_id = curr_otu_id
    
                for seq_id in otu_info[1:]:
                    sequence = self.retrieve_seq_by_id(seq_id,seqs)
                    if sequence is None:
                        logging.error('Unable to retrieve {}'.format(seq_id))
                        raise Exception('divvy: Unable to retrieve '+seq_id)
                    divvied_seqs += sequence
                    seqs_read += 1
                if seqs_read >= seqs_per_file:
                    logger.info("Writing at {} over {} sequences".
                        format(seqs_read,seqs_per_file))
                    out_filename = "{}_{}_{}{}".format(seqs_name, base_otu_id,
                        curr_otu_id, seq_ext)
                    out_filename = os.path.join(divvy_dir,out_filename)
                    # Write and reset all three values
                    logger.info("divvy: writing {}...".format(out_filename))
                    with open(out_filename,'w') as out_file:
                        out_file.write(divvied_seqs)
                    base_otu_id,divvied_seqs,seqs_read = None,"",0
                    file_count +=1 
            else: # write last file
                logger.info("Writing at {} over {} sequences".
                        format(seqs_read,seqs_per_file))
                out_filename = "{}_{}_{}{}".format(seqs_name, base_otu_id,
                        curr_otu_id, seq_ext)
                out_filename = os.path.join(divvy_dir,out_filename)
                # Write and reset all three values
                logger.info("divvy: writing {}...".format(out_filename))
                with open(out_filename,'w') as out_file:
                    out_file.write(divvied_seqs)
                base_otu_id,divvied_seqs,seqs_read = None,"",0
                file_count +=1 
    
                return file_count

    def sort_seqs_by_otu(self,otu_map,seq_file,seqs_per_file,divvy_output_dir='divvy'):
        """
        Scan the otu table, writing associated sequences to a new
        .fna file. Attempt to assign roughly equivalent numbers of 
        sequences per file, by taking the total number of sequences
        divided by the number of processes to be run.
    
        Files will be output to the output/divvy_seqs folder.
        Filenames use the sequence filename before the '.fna',
        with the range of OTU numbers suffixed on.
        For example, a divvied file with sequences from OTUs 1 to 122
        and a base sequence filename of seqs_R1_001.fna reads as follow:
        seqs_R1_001_1_122.fna
        """ 
        logger.info("divvying initiated")
        logger.debug("OTU map is "+otu_map)
        logger.debug("Sequence file located at "+seq_file)
        logger.debug("Target number of sequences per file: "+str(seqs_per_file))
        logger.debug("Target directory for divvied sequences: "+divvy_output_dir)
    
        # Check directory to store these sequences
        divvy_dir = self.validate_directory(divvy_output_dir)
        logger.debug('directory validated: '+divvy_dir)
    
        # divvy OTUs into separate files
        count = self.divvy_otus(divvy_dir,otu_map,seq_file,seqs_per_file)
    
        #return directory of sorted sequences
        logger.info("divvying complete")
        return divvy_dir,count
    

class Consolidate:
    def __init__(self):
        self.logger = logging.getLogger('Consolidate')

    def concat_otus(self,otu_dir,output_dir,seq_file):
        """
        Compiles the n OTU tables into one file and writes it 
        to hc_{seq filename}_otus.txt
        Returns a tuple of the consolidated OTU map and the number 
        of OTUs written.
        """
        logger.info("Consolidation initiated")
    
        # Build filepath
        logger.debug("seq_file is " + seq_file)
        seqs_path,seqs_filename = os.path.split(seq_file)
        logger.debug(seqs_path + ", " + seqs_filename)
        seqs_basename,seqs_ext = os.path.splitext(seqs_filename)
        logger.debug(seqs_basename + ", " + seqs_ext )
        otu_filename = 'hc_{seqs}_otus.txt'.format(seqs=seqs_basename)
        otu_dest = os.path.join(output_dir,otu_filename)
        logger.info("Consolidating output at " + otu_dest)
    
        with open(otu_dest,'w') as otu_file:
            otu_num = 0
            # get all otu filepaths     
            tmp_dir = os.path.join(otu_dir,"*/*.txt")
            dir_list = glob.glob(tmp_dir)
            logger.debug("Consolidating these files:\n" + str(dir_list))
            logger.debug("Consolidating files under " + otu_dir )
    
            for f_path in dir_list:
                with open(f_path,'r') as curr_f:
                    logger.info("Consolidating " + curr_f.name)
                    files_otus = 0
                    for line in curr_f:
                        otu_seqs = line.split('\t',1)[1]
                        to_write = str(otu_num)+' '+ otu_seqs
                        otu_file.write(to_write)
                        otu_num += 1
                        files_otus += 1
                    logging.debug('Contains {} OTUs'.format(files_otus))
        logger.info("Wrote {} OTUs to {}".format(otu_num,otu_dest))
        return otu_dest,otu_num 
    

class MPICluster:
    def __init__(self,script_name='mPyClust.py'):
        self.logger = logging.getLogger('enterMPY')
        self.py_script = script_name

    def mpi_cluster(self,out_dir,seq_dir, num_procs, sim):
        logger.info('Entering second clustering')
        logger.info('Divvied sequence directory given as ' + seq_dir)
        logger.info(str(num_procs) + " processes requested")
    
        # Create hc_{sim} parent directory
        parent_dir = os.path.join(out_dir,'hc_{}'.format(str(int(sim*100))))
        try:
            os.makedirs(parent_dir)
        except OSError as e:
            logger.warning('Overwrite will occur on parent directory\n\t'+ str(e))
    
        passed_args = ['mpiexec',
                       '-n', str(num_procs),
                       'python', self.py_script, parent_dir, seq_dir, str(sim)]
        logging.debug('Subprocess arguments: {}'.format(passed_args))
        subprocess.call(passed_args,shell=False)
    
        logging.info('Parent directory for second-cluster: {}'.format(parent_dir))
        return parent_dir


def validate_args():
    if not os.path.exists(args.sequences):
        logger.error("Sequence file does not exist")
        raise Exception("Sequence file does not exist")

def retrieve_otu_map_filename(otu_dir):
    otu_map = glob.glob(os.path.join(otu_dir,'*_otus.txt'))
    if otu_map is None:
        logger.error('OTU map filename could not be retrieved')
        raise Exception('Failed to locate otu map from first clustering')
    logger.info("OTU_map is {}".format(otu_map[0]))
    return otu_map[0]


def count_sequences(seq_file):
    #Get number of sequences in the file
    grep_args = ['grep','-c',">",seq_file]
    proc = subprocess.Popen(grep_args,stdout=subprocess.PIPE,shell=False)
    num = int(proc.stdout.read().rstrip('\n'))
    logger.info('{} sequences counted'.format(num))
    return num

def log_uclust_otu_count(uclust_dir):
    """
    Retrieves the OTU count from the .log file generated by uclust
    Returns the number created or -1 if it hits a problem
    """
    log_files = glob.glob(os.path.join(uclust_dir,'*.log'))
    with open(log_files[0],'r') as uclust_log:
        num_otus = -1
        for line in uclust_log:
            if "Num OTUs" in line:
                try:
                    num_otus = int(line.split(':')[1])
                    logger.info('{} OTUs created in first cluster'.format(num_otus))
                    break
                except Exception:
                    logger.warning('Number of OTUs could not be retrieved')
    return num_otus

def setup_dir(_top='output'):
    global num_seqs
    top_dir = _top
    # Build run-specific hiClust directory
    if num_seqs is 0:
        num_seqs = count_sequences(args.sequences)
    hc = 'hc_{}_{}_{}'.format(first_num,second_num,num_seqs)
    hc = os.path.join(top_dir,hc)
    divvy = os.path.join(hc,'divvy')
    logger.info('output directories:\n'+
        top_dir + '\n' + hc + '\n' + divvy )

    try: 
        if not os.path.exists(divvy):
            # makedirs builds needed intermediate directories
            os.makedirs(divvy)
    except OSError:
        logger.warning('error creating directories')

    return top_dir,hc,divvy

def main():
    validate_args()
    num_seqs = count_sequences(args.sequences)
    top_dir,hc_dir,divvy_dir = setup_dir(args.output_directory)
    uclust_dir = '{}_uclust_otus'.format(first_num)
    uclust_dir = os.path.join(hc_dir,uclust_dir)
    logger.info('uClust directory at '+ uclust_dir)

    # 1) Cluster at 90% 
    if not args.skipUclust:
        uclust_args = [method,
            '-i','{}'.format(args.sequences),
            '-o',uclust_dir,
            '-s',str(sim_first),
            '-A']
        subprocess.call(uclust_args,shell=False)    
    logger.info('First clustering complete')
    log_uclust_otu_count(uclust_dir)

    # 2) Separate sequences for re-clustering
    # Build OTU map filename from uclust directory
    otu_map = retrieve_otu_map_filename(uclust_dir)
    # Split sequences into files by first-cluster relationships
    divvier = Divvy()
    (seq_dir, num_files) = divvier.sort_seqs_by_otu(otu_map,
                                                    args.sequences,
                                                    num_seqs/args.processors,
                                                    divvy_dir)  
    logger.info("Divvying returned {} files at {}".format(num_files,seq_dir))

    # 3) Re-cluster in parallel
    # Ultimately we need to run as many processes as we have divvied 
    # sequence files - should be the same, however.
    mpiCluster = MPICluster()
    hc_otu_dir = mpiCluster.mpi_cluster(hc_dir,seq_dir,num_files,sim_second)
    logger.info('Second clustering complete')

    # 4) Consolidate sequences by 
    consolidator = Consolidate()
    out_out_dir,final_otu_count = consolidator.concat_otus(hc_otu_dir,hc_dir,args.sequences)
    logger.info('OTUs consolidated')
    logger.info(str (final_otu_count) + " OTUs created")

    print "OTUs",final_otu_count

    if not args.keep:
        rm_args = ['rm','-rf',uclust_dir,divvy_dir,hc_otu_dir]
        logger.info('Deleting intermediate files')
        subprocess.call(rm_args,shell=False)


if __name__ == '__main__':
    main()