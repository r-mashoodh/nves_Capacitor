QClus v0.1.0 installed from GitHub.

Steps to reproduce results:

- Run the computation of unspliced fraction per barcode:


	`mkdir QClusOut_NC_nves`

	in python:

	```python
	import pandas as pd
	import qclus as qc
	import scanpy as sc

	#TODO define cellranger_dir and outdir

	# run
	bam = f"{cellranger_dir}/possorted_genome_bam.bam"
	bam_index = f"{cellranger_dir}/possorted_genome_bam.bam.bai"
	barcodes_path = f"{cellranger_dir}/raw_feature_bc_matrix/barcodes.tsv.gz"

	fraction_unspliced = qc.utils.fraction_unspliced_from_bam(bam_path=bam, bam_index_path=bam_index, barcodes_path=barcodes_path, cores=24)
	fraction_unspliced.to_csv(f"{outdir}/fraction_unspliced.csv")
	```

- Blast to beetle the qclus reference set of human genes highly expressed in nucleus:

	`diamond makedb --threads 32 --in nves_longestT_pep.ok.fa --db nves_longestT_pep.dmnd`

	`diamond blastp --query human_prot_nucl_qclus.fasta --db nves_longestT_pep.dmnd --max-hsps 1 --evalue 1e-5 --outfmt 6 --out nves_pep.vs.human_nucl_prot_qclus.blp -p 36`

- Get gene name of nves genes that are best blast hit (e-val<10**-10):

	in python:

	```python
	with open('nves_pep.vs.human_nucl_prot_qclus.blp', 'r') as infile:
		res = {}
		for line in infile:
			human_gene, nves_gene = line.strip().split('\t')[:2]
			evalue = float(line.strip().split('\t')[-2])
			if evalue <= 10**-10 and human_gene not in res:
				res[human_gene] = nves_gene.split('.t')[0]


	with open('table_Nvi2_genes.tsv', 'r') as infile:
		my_nves_genes = []
		for line in infile:
			g, gname = line.strip().split('\t')[:2]
			if g in res.values():
				my_nves_genes.append(gname)

	print(my_nves_genes)
	```

	This prints the list below which to copy-paste into the qclus gene_lists.py (${CONDA_PREFIX}/lib/python3.11/site-packages/qclus/gene_lists.py) file: 
	['FOXP1', 'OSA', 'BNC2', 'FBXL7', 'ITPR', 'GEPH', 'RERE', 'TTC28', 'HR3', 'WWOX', 'TBC13', 'IMP2L', 'MED13', 'KGP25', 'EXOC6', 'EXD-3', 'EXOC4', 'SHEP', 'N42L2', 'EXT2', 'SYNE1-2', 'PARD3']

	Do the same for MT gene names (`grep '>MT-' nves_longestT_pep.fa`)

- Run qclus filtering:
	 `python filter_nuclei_qclus.py --unspliced 'QClusOut_NC_nves/fraction_unspliced.csv' --cellranger_dir '../expression_data/single-nucleus/cellranger/NC_nves/outs/' --use_matrix 'raw' --outdir 'QClusOut_NC_nves/' --min_genes 200`

	 `python filter_nuclei_qclus.py --unspliced 'QClusOut_FC_nves/fraction_unspliced.csv' --cellranger_dir '../expression_data/single-nucleus/cellranger//FC_nves/outs/' --use_matrix 'raw' --outdir 'QClusOut_FC_nves/' --min_genes 215`



