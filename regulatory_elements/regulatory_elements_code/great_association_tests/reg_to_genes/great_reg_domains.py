"""
Script for regions to gene association using the "GREAT" approach.

Command for beetle: python great_reg_domains.py --genes allgenes.wn.bed --chr genome_data/Ni2.chrsize.txt --output nicves_reg_domains.bed --format genes
"""

import sys

import argparse

from collections import OrderedDict

def load_chr(input_file):
    d = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            chrom, size = line.strip().split("\t")
            d[chrom] = int(size)
    return d


def load_genes(input_file, strip_name = False, gformat='tss'):
    d = OrderedDict()
    with open(input_file, 'r') as infile:
        for line in infile:
            if gformat == 'tss':
                chrom, tss, _, gene, _, strand = line.strip().split("\t")[:6]
                tss = int(tss)
            else:
                chrom, start, end, gene, _, strand = line.strip().split("\t")[:6]
                assert strand in ["1", "-1", '-', '+']
                if strand == '1' or strand == '+':
                    tss = int(start)
                else:
                    tss = int(end) - 1
            if strip_name:
                gene = '.'.join(gene.split('.')[:-1])
            d[gene] = chrom, tss, strand
    d = dict(sorted(d.items(), key=lambda item: (item[1][0], item[1][1])))
    return d


def basal_domains(genes, chr_sizes, up=5000, down=1000): #5 kb upstream #1kb downstream
    domains = {}
    for gene in genes:
        chrom, tss, strand = genes[gene]
        assert strand in ["1", "-1", '-', '+']
        if strand == '1' or strand == '+':
            domain_start = max(tss-up, 0) 
            domain_end = min(tss+down, chr_sizes[chrom]) 
        else:
            domain_end = min(tss+up, chr_sizes[chrom])
            domain_start = max(tss-down, 0)

        domains[chrom] = domains.get(chrom, OrderedDict())
        domains[chrom][(domain_start, domain_end)] = domains[chrom].get((domain_start, domain_end), [])
        domains[chrom][(domain_start, domain_end)].append(gene)
    for chrom in domains:
        tmp = domains[chrom]
        tmp = OrderedDict(sorted(tmp.items(), key=lambda item: (item[0][0], item[0][1])))
        domains[chrom] = tmp
    return domains


def basal_extension(domains, chr_sizes, max_extension=1000000):
    extended_domains = {}
    for chrom in domains.keys():
        extended_domains[chrom] = {}
        intervals = list(domains[chrom].keys())
        for i, interval in enumerate(intervals):
            start, end = interval

            # print(chrom, interval)

            next_i = i + 1
            prev_i = i - 1

            if prev_i > 0:
                prev_end = intervals[prev_i][1]
                if prev_end < start:
                    start = max(prev_end, start - max_extension) #max extension 1Mb


            else:
                start = max(start - max_extension, 0) #max extension 1Mb

            if next_i < len(intervals):
                next_start = intervals[next_i][0]
                if next_start > end:
                    end = min(end + max_extension, next_start) #max extension 1Mb

            else:
                end = min(end + max_extension, chr_sizes[chrom]) #max extension 1Mb
            genes = domains[chrom][interval]
            extended_domains[chrom][(start, end)] = extended_domains.get((start, end), [])
            # print(start, end)
            for gene in genes:
                extended_domains[chrom][(start, end)].append(gene)

                # if i<10:
                #     print(gene, chrom)
                #     print((start, end))
    return extended_domains

            

def write_out(domains, output_file):
    with open(output_file, 'w') as out:
        for chrom in domains:
            for dom in domains[chrom]:
                for gene in domains[chrom][dom]:
                    out.write('\t'.join([chrom, str(dom[0]), str(dom[1]), gene])+'\n')

if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('--genes', type=str, required=True, 
                        help="5-columns .bed file with genes TSS (chr, tss, tss+1, gene, strand), WARNING:: the code assumes the bed is sorted !!!")

    PARSER.add_argument('--chr', type=str, required=True, help="file with 2 col: chr_name and chr_size")

    PARSER.add_argument('--output', type=str, required=True)

    PARSER.add_argument('--up', type=int, required=False, default=5000)

    PARSER.add_argument('--down', type=int, required=False, default=1000)

    PARSER.add_argument('--max_extension', type=int, required=False, default=1000000)

    PARSER.add_argument('--strip', required=False, action="store_true")

    PARSER.add_argument('--no_extension', required=False, action="store_true")

    PARSER.add_argument('--format', required=False, default='tss') #can be one of tss or genes

    ARGS = vars(PARSER.parse_args())

    if ARGS['up'] == 5000 and ARGS['down'] == 1000:
        sys.stderr.write("Assigning regulatory domains to genes, using GREAT default association: basal 5kb upstream, 1kb downstream + extension.\n")
    else:
        u = ARGS['up'] 
        d = ARGS['down'] 
        sys.stderr.write(f"Assigning regulatory domains to genes, using custom GREAT association: basal {u/1000}kb upstream, {d/1000}kb downstream + extension.\n")

    sys.stderr.write("WARNING: the script assumes the input bed with tss coord is sorted, please run bedtools sort first. (Ignore the warning if you did).\n")

    assert ARGS['format'] in ['genes', 'tss']

    CHR_SIZE = load_chr(ARGS["chr"])
    GENES = load_genes(ARGS["genes"], ARGS['strip'], gformat=ARGS['format'])

    BASAL = basal_domains(GENES, CHR_SIZE, up=ARGS['up'], down=ARGS['down'])

    if not ARGS['no_extension']:
        REG_DOMAINS = basal_extension(BASAL, CHR_SIZE, max_extension=ARGS['max_extension'])
    else:
        REG_DOMAINS = BASAL

    write_out(REG_DOMAINS, ARGS["output"])

    sys.stderr.write(f"DONE! output saved in {ARGS['output']}.\n")


