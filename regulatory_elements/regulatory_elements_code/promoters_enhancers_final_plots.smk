'''
Snakemake pipeline for the defintion of N. vespilloides promoters and enhancers using ATAC and H3K4me3 enrichment:

- definition of enhancers and promoters
- distance to tss plots
- differential accessibility plots (same commands for differential binding or activity)

Was run using the '../envs/snake.yaml environment'
'''

ATAC_peaks = 'atac_peaks_Genrich_all.bed'
H3k4me3 = '../peaks/H3K4me3.q18.narrowPeak'
GTF = '../../genome_assembly/genome_files/genes/cellranger_gtf/Nvi2.liftoff.mitohifi.final.ok.genename.cellranger.gtf'

OUT = 'output'

BAM_DIR = 'bams/'

SAMPLES, = glob_wildcards(BAM_DIR+'{samples}.sorted.bam')
FC = sorted([i for i in SAMPLES if 'FC' in i])
NC = sorted([i for i in SAMPLES if 'NC' in i])
SAMPLES = FC + NC

ATAC = ['A3_FC',
'A4_FC',
'A5_FC',
'A6_FC',
'A3_NC',
'A4_NC',
'A5_NC',
'A6_NC']


rule all:
	input: OUT + "/plots/heatmaps/cat_all.DAreg.log.svg",
	       OUT + "/plots/tss_dist.svg",
	       OUT + "/plots/heatmaps/delta_vector.txt"


rule intersect:
	'''
	Define promoters by overlap of H3K4me3 peaks and atac peaks.
	'''
	input: a = ATAC_peaks, b = H3k4me3
	output: OUT + '/promoters.bed'
	shell: 'bedtools intersect -a {input.a} -b {input.b} -wa | uniq > {output}'

rule enhancers:
	'''
	Putative enhancers = atac peaks without H3K4me3.
	'''
	input: p = OUT + '/promoters.bed', r = ATAC_peaks
	output: e = OUT + '/enhancers.bed', p = OUT + '/promoters_ids.txt'
	shell: 'cut -f 4 {input.p} > {output.p} && grep -Fwvf {output.p} {input.r} > {output.e}'

rule get_tss_from_gtf:
	'''
	Extract TSS location from the gtf
	'''
	input: gtf = GTF
	output: OUT + "/TSS/tss_all_location.bed"
	shell: "python get_tss_from_gtf.py -i {input.gtf} -o {output}_temp.bed --allTSS && bedtools sort -i {output}_temp.bed > {output} && rm {output}_temp.bed"


rule plot_elements:
	'''
	QC plots for enhancers and promoters
	'''
    input:
        infiles = expand(OUT + '/{reg}.bed', reg=['promoters', 'enhancers']),
        tss = OUT + "/TSS/tss_all_location.bed"
    output:
        bar = OUT + "/plots/barplot_nb_reg.svg",
        bx =  OUT + "/plots/boxplot_size_reg.svg",
        tss = OUT + "/plots/tss_dist.svg"  

    shell:
        "python control_plots.py -i {input.infiles} -l promoters enhancers "
        "-t 'beetle' -o {output.bar} {output.bx} {output.tss} "
        "-tss {input.tss}"


rule filter_bam_strict:
	'''
	filter bam for mapq>13
	'''
	input: "bams/{sample}.sorted.bam"
	output: OUT + "/filtered_bams/{sample}_mapq13.bam"
	conda: '../envs/samtools.yaml'
	shell: "samtools view -f 2 -q 13 {input} -o {output} && samtools index {output}"


rule coverage:
	'''
	bam to bigwig
	'''
	input: b = OUT + "/filtered_bams/{sample}_mapq13.bam"
	output: OUT + '/bw/{sample}.bw'
	threads: 16
	conda: '../envs/deeptools.yaml'
	shell: "bamCoverage -b {input.b} -o {output} -p {threads} --normalizeUsing RPKM"

rule element_sets:
	'''
	Subclassify differential peaks into differential promoters and differential enhancers
	'''
	input:
		bed1 = OUT + '/{reg}.bed',
		bed2 = '{UpCat}.bed'
	output: OUT + '/{UpCat}-{reg}.bed'
	shell: 'bedtools intersect -a {input.bed2} -b {input.bed1} -wa > {output}'

rule bw_to_matrix:
	'''
	Deeptools matrix for heatmap plot
	'''
	input: 
		bed = expand(OUT + '/{UpCat}-{reg}.bed', UpCat=['DA_UpCare', 'DA_UpNoCare'], reg=['promoters', 'enhancers']),
		bw = expand('bw/avg_atac_{c}.bw', c=['FC', 'NC'])#bw = expand('../atac/out_beetle_atac/bw/{atac}.bw', atac=ATAC) #expand(OUT + '/bw/{sample}.bw', sample=SAMPLES) + 
	output: OUT + "/plots/heatmaps/cat_all.DAreg.matrix.gz"
	conda: '../envs/deeptools.yaml'
	threads: 16
	shell:
		"computeMatrix reference-point --referencePoint center -R {input.bed} "
		"-S {input.bw} -o {output} --beforeRegionStartLength 1000 --afterRegionStartLength 1000 -p {threads}"

# rule avg_fc:
# 	input:
# 		chrsize = 'Ni2.chrsize.txt',
# 		bw = expand('../atac/out_beetle_atac/bw/{atac}.bw', atac=[i for i in ATAC if 'FC' in i])
# 	output:
# 		temp(OUT + '/bw/avg_atac_FC.bedGraph')
# 	shell: "wiggletools write_bg {output} mean {input.bw}"

# rule avg_nc:
# 	input:
# 		chrsize = 'Ni2.chrsize.txt',
# 		bw = expand('../atac/out_beetle_atac/bw/{atac}.bw', atac=[i for i in ATAC if 'NC' in i])
# 	output:
# 		temp(OUT + '/bw/avg_atac_NC.bedGraph')
# 	shell: "wiggletools write_bg {output} mean {input.bw}"

# rule tobw:
# 	input: chrsize = 'Ni2.chrsize.txt', bg = OUT + '/bw/avg_atac_{c}.bedGraph'
# 	output:  bw = OUT + '/bw/avg_atac_{c}.bw', tmp = temp(OUT + '/bw/avg_atac_{c}.ok.bedGraph')
# 	shell: 'grep -v chrMT {input.bg} > {output.tmp} && bedGraphToBigWig {output.tmp} {input.chrsize} {output.bw}'


rule delta:
	'''
	Create delta track
	'''
	input: OUT + "/plots/heatmaps/cat_all.DAreg.matrix.gz"
	output: OUT + "/plots/heatmaps/delta_vector.txt"
	run:
		import gzip
		import json
		import numpy as np
		import pandas as pd


		def get_center(matrix):
			num_bins = matrix.shape[1]
			mid = num_bins // 2
			if num_bins % 2 == 0:
				# Even: Average the 10 central bins (5 left, 5 right of center)
				return np.mean(matrix[:, mid-25 : mid+25], axis=1) #matrix[:, mid-5 : mid+5]
			else:
				# Odd: Average the 9 central bins (center bin + 4 left + 4 right)
				return np.mean(matrix[:, mid-24 : mid+25], axis=1)


		with gzip.open(input[0], 'rt') as f:
			header_line = f.readline()
			header = json.loads(header_line[1:])
			bounds = header['sample_boundaries'] # [0, 200, 400]

			data = pd.read_csv(f, sep="\t", header=None)
			
			# Load numerical data (skipping first 6 metadata columns)
			# data = np.loadtxt(f, delimiter='\t')

			signal = np.array(data.iloc[:, 6:])

			vals_fc = signal[:, bounds[0]:bounds[1]]
			vals_nc = signal[:, bounds[1]:bounds[2]]

			vec_a = get_center(vals_fc)
			vec_b = get_center(vals_nc)

			# Delta: (A-B) / max(A,B)
			denom = np.maximum(vec_a, vec_b)
			# delta = (vec_a - vec_b)
			delta = np.divide(vec_a - vec_b, denom, out=np.zeros_like(vec_a))

			np.savetxt(output[0], delta, fmt='%g')



rule log2:
	'''
	log2 rpkm
	'''
	input: OUT + "/plots/heatmaps/cat_all.DAreg.matrix.gz"
	output: OUT + "/plots/heatmaps/cat_all.DAreg.log.matrix.gz"
	run:
		import gzip
		import pandas as pd
		import numpy as np

		with gzip.open(input[0], 'rb') as f:
		    line = f.readline()
		    header = line.decode('utf-8')
		    data = pd.read_csv(f, sep="\t", header=None)

		data.iloc[:, 6:] = np.log2(data.iloc[:, 6:] + 1)

		with gzip.open(output[0], 'wb') as f:
		    f.write(header.encode('utf-8'))
		    data.to_csv(f, sep="\t", header=False, index=False)

rule plot_heatmap:
	'''
	Hetmap plot
	'''
	input: OUT + "/plots/heatmaps/cat_all.DAreg.log.matrix.gz"
	output: OUT + "/plots/heatmaps/cat_all.DAreg.log.svg"
	conda: '../envs/deeptools.yaml'
	shell: "plotHeatmap -m {input} -o {output} --regionsLabel 'Care-prom' 'Care-enh' 'NoCare-prom' 'NoCare-enh' --zMin 3 --zMax 8 --plotType se "
	       "--refPointLabel 0 --colorMap BuPu --interpolationMethod nearest --heatmapHeight 12 --heatmapWidth 3" #yMax 11 for enh 



