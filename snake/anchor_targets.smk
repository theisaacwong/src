import re
import pandas as pd
import json

configfile: "target_sites.yaml"

MANIFEST = config.get("NEWHAPS")
BED = config.get("BED")
REF = config['REF']

GENE = 'RPTOR'

bed_df = pd.read_csv(BED, sep="\t", header=None, usecols=[0,1,2,3], names=["chr", "start", "end", 'NAME'], dtype=str)
bed_df = bed_df.loc[bed_df['NAME'] == GENE]
bed_df.set_index("NAME", inplace=True)

anchor_df = pd.read_csv(f'results/{GENE}/anchor_manifest.tab', sep='\t', index_col=['ANCHOR_ID'])
manifest_df = pd.read_csv(MANIFEST, sep="\t", index_col=["SAMPLE"])


def find_anchor(wildcards):
    return anchor_df.at[wildcards.anchor, 'FA']

def find_asm(wildcards):
    return manifest_df.at[wildcards.sample, 'HAP']

def cigar_tuple(cigar):
    """
    Takes a cigar string as input and returns a cigar tuple
    """
    opVals = re.findall(r'(\d+)([\w=])', cigar)
    lengths = [int(opVals[i][0]) for i in range(0, len(opVals))]
    ops = [opVals[i][1] for i in range(0, len(opVals))]
    tuples = []
    for i in range(len(lengths)):
        tuples.append([int(lengths[i]), ops[i]])
    return tuples

rule all:
    input:
        # expand("results/{region}/anchors/{sample}_{anchor}_asm_aln.paf", sample=manifest_df.index, region=bed_df.index, anchor=anchor_df.index),
        expand("results/{region}/asm_anchors_stats.tsv", region=bed_df.index)


rule anchor_aln:
    input:
        asm=find_asm,
        anchor=find_anchor
    output:
        paf="results/{region}/anchors/{sample}_{anchor}_asm_aln.paf"
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.12.0",
        "minimap2/2.28"
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 6
    shell:
        """
        minimap2 -x asm20 -c -t {threads} --eqx {input.asm} {input.anchor} > {output.paf}
        """

# rule anchor_local_aln:
#     input:
#         asm='results/{region}/fa/{sample}.fa',
#         anchor=find_anchor
#     output:
#         paf="results/{region}/anchors/{sample}_{anchor}_aln.paf"
#     envmodules:
#         "modules",
#         "modules-init",
#         "modules-gs/prod",
#         "modules-eichler/prod",
#         "miniconda/4.12.0",
#         "minimap2/2.28"
#     resources:
#         mem=10,
#         hrs=24,
#         disk_free=1,
#     threads: 6
#     shell:
#         """
#         minimap2 -x asm20 -c -t {threads} --eqx {input.asm} {input.anchor} > {output.paf}
#         """

rule summarize_asm_aln_anchors:
    input:
        pafs = expand("results/{region}/anchors/{sample}_{anchor}_asm_aln.paf",
            sample=manifest_df.index,anchor=anchor_df.index,region=bed_df.index),
    output:
        tsv = "results/{region}/asm_anchors_stats.tsv",
        range = "results/{region}/asm_anchors_range_variance.tsv"
    resources:
        mem=24,
        hrs=2,
    run:
        fout = open(output.tsv, 'w')
        anchor_dict = {}
        contig_sample_map = {}
        for paf in input.pafs:
            sample = paf.split('/')[-1].split('_')[0]
            for paf_line in open(paf):
                line = paf_line.strip().split()

                qchr = line[0]  ## CHM13 query
                anchor_id = qchr
                rchr = line[5]  ## other genomes
                contig_sample_map[rchr] = sample
                qstart = int(line[1])
                qend = int(line[2])

                rstart = int(line[7])
                rend = int(line[8])
                rlen = abs(rend - rstart) + 1
                qlen = abs(qend - qstart) + 1

                cg = [i.split(":")[-1] for i in line[12:] if i[:2] == 'cg'][0]

                cg_tuple = cigar_tuple(cg)
                iden = round((sum([int(i[0]) for i in cg_tuple if i[1] == '=' or i[1] == 'M']) / sum(
                    [int(i[0]) for i in cg_tuple if i[1] in {'=', 'M', 'X', 'D', 'I'}])) * 100,2)
                qCov = round(float(int(qend) - int(qstart)) / int(qlen) * 100,4)
                print(f'{rchr}\t{rstart}\t{rend}\t{iden}\t{qCov}\t{anchor_id}',file=fout)

                if rchr not in anchor_dict:
                    anchor_dict[rchr] = []
                anchor_dict[rchr].append([rstart, rend, iden])

        popdict = json.load(open('/net/eichler/vol28/projects/medical_reference/nobackups/Samples/igsr_samples_pop.json'))
        frange = open(output.range, 'w')
        print('Contig\tLength\tAnchor1_start\tAnchor1_end\tAnchor1_iden\tAnchor2_start\tAnchor2_end\tAnchor2_iden\tPop\tSample', file=frange)
        for rchr, anchor_list in anchor_dict.items():
            if len(anchor_list) == 2:
                sample = contig_sample_map[rchr]
                pop = popdict[sample] if sample not in ['PTR', 'PPA', 'GGO', 'PPY', 'PAB', 'SSY'] else sample
                sorted_anchors = sorted(anchor_list, key=lambda x:x[0])
                length = sorted_anchors[1][0] - sorted_anchors[0][1]
                anchor1 = '\t'.join([str(ele) for ele in sorted_anchors[0]])
                anchor2 = '\t'.join([str(ele) for ele in sorted_anchors[1]])
                print(f'{rchr}\t{length}\t{anchor1}\t{anchor2}\t{pop}\t{sample}', file=frange)

