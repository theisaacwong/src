#new version for trio variant calling
#October 2021

import pysam
import pandas as pd
import argparse
from tqdm import tqdm
from joblib import Parallel, delayed

M=0 #M  BAM_CMATCH      0
I=1 #I  BAM_CINS        1
D=2 #D  BAM_CDEL        2
N=3 #N  BAM_CREF_SKIP   3
S=4 #S  BAM_CSOFT_CLIP  4
H=5 #H  BAM_CHARD_CLIP  5
P=6 #P  BAM_CPAD        6
E=7 #=  BAM_CEQUAL      7
X=8 #X  BAM_CDIFF       8
B=9 #B  BAM_CBACK       9
NM=10 #NM       NM tag  10
conRef  =       [M, D, N, E, X] # these ones "consume" the reference
conQuery=       [M, I, S, E, X] # these ones "consume" the query
conAln  =       [M, I, D, N, S, E, X] # these ones "consume" the alignments

def get_bam(sample, family):
	bam_path = bams.loc[(bams['sample'] == sample) & (bams['family'] == family), 'cram'].values[0]
	bam = pysam.AlignmentFile(bam_path)
	return bam

def assign_base(allele):
	if allele == "A":
		base = 0
	elif allele == "C":
		base = 1
	elif allele == "G":
		base = 2
	elif allele =="T":
		base = 3
	else:
		base = 4

	return(base)

def find_coverage(chr, pos, bam):
	#use pysam to count base calls at site, then format as list
	cov = bam.count_coverage(chr, start=pos-1, stop=pos, quality_threshold = None)
	cov = [  x[0] for x in cov ]

	total_cov = 0

	for pileupcolumn in bam.pileup(chr, start=pos-1, stop=pos):
		if pileupcolumn.pos == pos:
			total_cov = pileupcolumn.n

	return cov, total_cov

def base_counts(ref_base, alt_base, sample_cov):
	ref_count = sample_cov[ref_base]
	alt_count = sample_cov[alt_base]

	counts = str(ref_count) + "," + str(alt_count)

	return counts

def indel_count(bam, chr, pos, alt, indel_length):
	alt_count = 0
	num_reads = 0

	for read in bam.fetch(chr, pos-1, pos):
		num_reads = num_reads + 1
		ref_idx = pos - read.reference_start
		ref_pos = 0
		read_pos = 0

		cigartuples = read.cigartuples
		if cigartuples is None:
			continue
		for op, length in cigartuples:
			if op in conRef:
				ref_pos = ref_pos + length
			if op in conQuery:
				read_pos = read_pos + length
			#if deletion, it will be indexed at the end of the deleted sequence
			if ref_pos == ref_idx+length and op == 2 and length == indel_length:
				if read.query_sequence[read_pos-1:read_pos] == alt:
					alt_count = alt_count + 1
			#if insertion, it will be indexed at the start of the inserted sequence
			if ref_pos == ref_idx and op == 1 and length == indel_length:
				if read.query_sequence[read_pos-length-1:read_pos] == alt:
					alt_count = alt_count + 1

	ref_count = num_reads - alt_count
	counts = str(ref_count) + "," + str(alt_count)

	return counts, num_reads

def fix_counts(row, sample, sample_count, sample_total):
	if sample == 'fa':
		idx = 9
	elif sample == 'mo':
		idx = 10
	else:
		idx = 11

	info = row[idx].split(':')
	info[1] = sample_count
	info[2] = str(sample_total)
	new_info = ':'.join(info)
	return new_info

def snv_depth(line, mom, dad, kid, sample_id):
	chr = line[0]
	pos = int(line[1])
	ref = line[3]
	alt = line[4]

	ref_base = assign_base(ref)
	alt_base = assign_base(alt)

	dad_cov, dad_total = find_coverage(chr, pos, dad)
	mom_cov, mom_total = find_coverage(chr, pos, mom)
	kid_cov, kid_total = find_coverage(chr, pos, kid)

	dad_count = base_counts(ref_base, alt_base, dad_cov)
	mom_count = base_counts(ref_base, alt_base, mom_cov)
	kid_count = base_counts(ref_base, alt_base, kid_cov)

	new_dad = fix_counts(line, 'fa', dad_count, dad_total)
	new_mom = fix_counts(line, 'mo', mom_count, mom_total)
	new_kid = fix_counts(line, sample_id, kid_count, kid_total)

	line[9:12] = new_dad, new_mom, new_kid

	return line

def indel_depth(line, mom, dad, kid, sample_id):
	chr = line[0]
	pos = int(line[1])
	ref = line[3]
	alt = line[4]
	indel_length = abs(len(ref) - len(alt))

	dad_count, dad_total = indel_count(dad, chr, pos, alt, indel_length)
	mom_count, mom_total = indel_count(mom, chr, pos, alt, indel_length)
	kid_count, kid_total = indel_count(kid, chr, pos, alt, indel_length)

	new_dad = fix_counts(line, 'fa', dad_count, dad_total)
	new_mom = fix_counts(line, 'mo', mom_count, mom_total)
	new_kid = fix_counts(line, sample_id, kid_count, kid_total)

	line[9:12] = new_dad, new_mom, new_kid

	return line


if __name__ == "__main__":

	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input", required=True, help="Candidate de novo file")
	ap.add_argument("-b", "--bams", required=True, help="Manifest of bam files")
	ap.add_argument("-s", "--sample", required=True, help="ID of sample for de novo calling")
	ap.add_argument("-f", "--family", required=True, help="ID of family")
	args = vars(ap.parse_args())

	bams = pd.read_csv(args["bams"], sep='\t')
	bams['family'] = bams['family'].astype('str')

	mom_bam = get_bam("mo", args['family'])
	dad_bam = get_bam("fa", args['family'])
	kid_bam = get_bam(args["sample"], args['family'])

	output =  open("%s.%s.depth.adjusted.candidate.denovo.txt"%(args["family"], args["sample"]), "w")

	with open(args['input']) as candidates:
		for line in candidates:
			line = line.split('\t')
			if line[4] == '*':
				continue
			if len(line[3]) == len(line[4]):
				new_line = snv_depth(line, mom_bam, dad_bam, kid_bam, args['sample'])
			else:
				new_line = indel_depth(line, mom_bam, dad_bam, kid_bam, args['sample'])

			output.write("\t".join(new_line))

	output.close()


