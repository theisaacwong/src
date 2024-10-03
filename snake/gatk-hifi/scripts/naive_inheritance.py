import sys
import gzip
import pandas as pd
import argparse
from tqdm import tqdm

def read_vcf(file):
    """
    Reads in gzipped vcf.
    :param file: Gzipped vcf file.
    :returns: List of header rows and list of variant rows.
    """
    header = []
    vcf = []

    with gzip.open(file, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                header.append(line)
            else:
                vcf.append(line.rstrip().split('\t'))
    return header, vcf


def assign_variant(row):
    """
    Assigns inheritance and transmission based on genotypes.
    :param row: Row of dataframe corresponding to variant.
    :returns: Transmission status ("TRANSMITTED=no/yes") and inheritance ("INH=mother_only/denovo_pro/etc."")
    """
    dad = row['fa'].split(":")[0].replace("|", "/")
    mom = row['mo'].split(":")[0].replace("|", "/")
    kid = row['kid'].split(":")[0].replace("|", "/")

    if dad == '0/0' and mom == '0/0' and kid == '0/0':
        return "TRANSMITTED=no", "INH=not_variant"
    elif dad == '0/0' and (mom == '0/1' or mom == '1/1') and kid == '0/0':
        return "TRANSMITTED=no", "INH=mother_only"
    elif (dad == '0/1' or dad == '1/1') and mom == '0/0' and kid == '0/0':
        return "TRANSMITTED=no", "INH=father_only"
    elif dad == '0/0' and mom == '0/0' and (kid == '0/1' or kid == '1/1'):
        return "TRANSMITTED=no", "INH=denovo_kid"
    elif (dad == '0/1' or dad == '1/1') and mom == '0/0' and kid == '0/1':
        return "TRANSMITTED=yes", "INH=fa_to_kid"
    elif dad == '0/0' and (mom == '0/1' or mom == '1/1') and kid == '0/1':
        return "TRANSMITTED=yes", "INH=mo_to_kid"
    elif dad == '1/1' and mom == '1/1' and kid == '1/1':
        return "TRANSMITTED=yes", "INH=unknown"
    elif (dad == '0/1' or dad == '1/1') and (mom == '0/1' or mom == '1/1') and (kid == '0/1' or kid == '1/1'):
        return "TRANSMITTED=yes", "INH=unknown"
    elif (dad == '0/1' or dad == '1/1') and (mom == '0/1' or mom == '1/1') and kid == '0/0':
        return "TRANSMITTED=no", "INH=parents_only"
    else:
        return "TRANSMITTED=unknown", "INH=unknown"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="vcf of all family variant calls")
    ap.add_argument("-c", "--caller", required=True, help="caller used to make variant calls")
    ap.add_argument("-f", "--family", required=True, help="family ID")
    args = vars(ap.parse_args())

    print("\nReading VCF...")

    header, vcf = read_vcf(args["input"])

    column_names = header[-1].rstrip().split('\t')
    df = pd.DataFrame(vcf)
    df.columns = column_names
    df.columns = df.columns.str.replace('%s.'%args["family"], '')

    sample = df.columns[-1]
    df.columns = [*df.columns[:-1], 'kid']

    transmission = []
    inheritance = []

    print("\nAssigning inheritance...")

    for index, row in tqdm(df.iterrows(), total = len(df)):
        tra, inh = assign_variant(row)
        row["INFO"] = row["INFO"] + ";" + tra + ";" + inh

    denovo = df[df["INFO"].str.contains("denovo_kid")]

    print("\nWriting output file...")

    caller = args["caller"]
    denovo.to_csv("%s.%s.candidate.denovo.sites.%s.txt"%(args["family"], sample,caller), sep = '\t', index = None, header = False)
