#/usr/bin/env/python

import argparse


def tss_longest_transcript(input_file, chlab='', coding_only=False, allTSS=False):
	d = {}
	k = 0
	with open(input_file, 'r') as infile:
		for line in infile:
			if line[0] == '#' or not line.strip():
				continue

			ch, _, feat, st, e, _, strand = line.split('\t')[:7]

			if chlab and chlab not in ch:
				continue

			if feat not in ['transcript', 'mRNA', 'tRNA', 'misc_RNA', 'ncRNA', 'rRNA']:
				continue

			if coding_only and 'protein_coding' not in line:
				continue

			assert strand in ['+', '-'], f"Error: strand should be + or -, {strand} found."

			gene_id = line.split(';')[0].split()[-1].strip('"')

			if strand == '+':
				bed_entry = '\t'.join([ch, str(int(st)-1), st, gene_id, '.', '+'])
			else:
				bed_entry = '\t'.join([ch, str(int(e)-1), e, gene_id, '.', '-'])

			if allTSS:
				d['tss'+str(k)] = ('NA', bed_entry)
				k += 1

			else:
				lg = int(e) - int(st) #may need +1 but it does not matter, we just use this to compare across isoforms
				if gene_id not in d or lg > d[gene_id][0]:

					d[gene_id] = (lg, bed_entry)

	return d


def write_bed(d, output_file):
	with open(output_file, 'w') as out:
		for gene in d:
			_, entry = d[gene]
			out.write(entry+"\n")



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input', required=True)

    parser.add_argument('-o', "--output", required=True)

    parser.add_argument('-c', "--chlab", required=False, default='')

    parser.add_argument("--coding_only", action='store_true')

    parser.add_argument("--allTSS", action='store_true') #get all tss not just longest transcript

    args = vars(parser.parse_args())

    tss = tss_longest_transcript(args['input'], args['chlab'], args['coding_only'], args['allTSS'])

    write_bed(tss, args['output'])