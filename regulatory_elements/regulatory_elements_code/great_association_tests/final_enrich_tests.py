import argparse
import os

import pandas as pd
import scipy.stats as ss


def hypergeom_genome(genes, regions, Ng):
    #number of genes in list
    Kg = len(genes)

    #number of genes associated with foreground regions
    ng = len({item for val in regions.values() for item in val})
    # print({item for val in regions.values() for item in val})

    #number of genes in gene list associated with foreground regions
    kg = len({item for val in regions.values() for item in val if item in genes})

    # len({regions[i] for i in regions if regions[i].intersection(genes)})
    # print(ng)
    # print('Overlap =', kg)
    # print({item for val in regions.values() for item in val if item in genes})


    pval = ss.hypergeom.sf(kg-1, Ng, Kg, ng)
    # print('Hypergeometric test p-value: ', pval)

    # print('Enrichment ratio: ', (kg/ng)/(Kg/Ng))
    # if pval < 0.05:
    #     print('#*#')
    return ng, kg


def binom_background(genes, regions, allpeaks):
    allpeaks['lg'] = allpeaks[2] - allpeaks[1]

    #total size of background regions with at least one gene
    tmp = allpeaks.drop_duplicates(3, keep='first')
    N = sum(tmp['lg'])

    ###total size of background associated with DEG
    tmp = allpeaks[allpeaks[5].isin(genes)]
    tmp = tmp.drop_duplicates(3, keep='first')
    K = sum(tmp['lg'])


    ##total size foreground regions: 
    tmp = allpeaks[allpeaks[3].isin(regions.keys())]
    tmp = tmp.drop_duplicates(3, keep='first')
    n = sum(tmp['lg'])

    #total size of foreground regions associated with DEG
    tmp = allpeaks[allpeaks[3].isin({i for i in regions if regions[i].intersection(genes)})]
    tmp = tmp.drop_duplicates(3, keep='first')
    k = sum(tmp['lg'])

    p = K/N

    # print('Fraction of genome associated with genelist: ', p)

    res = ss.binomtest(k, n, p)
    # print(k/n)
    # print('Bionomial test p-value: ', res.pvalue)
    # print('Enrichment ratio: ', k/n/p)

    # print(k/n)


def hyper_background(genes, allpeaks, ng, kg):
    
    Ng = len(allpeaks[5].unique())
    
    tmp = allpeaks[allpeaks[5].isin(genes)]
    Kg = len(tmp[5].unique())

    pval = ss.hypergeom.sf(kg-1, Ng, Kg, ng)
    print('Hypergeometric test p-value: ', pval)

    print('Enrichment ratio: ', (kg/ng)/(Kg/Ng))
    # print(pval)
    # if pval < 0.05:
    #     print('#*#')




if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    PARSER.add_argument('--genes', nargs='+', required=False, default=['../../../gene_expression/single_nucleus/DEG/Consistent_DEG_UpCare.txt', '../../../gene_expression/single_nucleus/DEG/Consistent_DEG_UpNoCare.txt']) #'Bulk_DEG_NoCareUP.txt', 'Bulk_DEG_CareUP.txt'

    PARSER.add_argument('--regions', nargs='+', required=False, default=['data/DHis_UpCareGREAT_ATACpeaksFINAL.enh.bed', 'data/DHis_UpNoCareGREAT_ATACpeaksFINAL.enh.bed']) #DHis_UpNoCareGREAT_ATACpeaksFINAL.bed 'DHis_UpCareGREAT_ATACpeaksFINAL.enh.bed', 'DHis_UpNoCareGREAT_ATACpeaksFINAL.enh.bed' DB_UPCare_to_genesGREAT.bed', 'DB_UPNoCare_to_genesGREAT.bed']) #['DA_UPCare_to_genesGREAT.bed', 'DA_UPNoCare_to_genesGREAT.bed']

    PARSER.add_argument('--allpeaks', type=str, required=False, default='data/H3K4me3_to_genes_GREAT_ATAC.enh.bed') #H3K4me3_to_genes_GREAT_ATAC.bed Regions_to_genes_sorted_GREAT.bed

    PARSER.add_argument('--regtype', type=str, required=False, default=None)

    PARSER.add_argument('--totalCodingGenes', type=int, required=False, default=12836)

    ARGS = vars(PARSER.parse_args())

    allpeaks = ARGS['allpeaks']
    allpeaks = pd.read_csv(allpeaks, sep='\t', header=None)
    if ARGS['regtype']:
        allpeaks = allpeaks[allpeaks[4] == ARGS['regtype']]

    for i in range(len(ARGS['genes'])):
        for j in range(len(ARGS['regions'])):
            genes = ARGS['genes'][i]
            regions = ARGS['regions'][j]


            print('***************************************************************')
            print(os.path.basename(genes), os.path.basename(regions))
            
            
            #pandas load and transform to list
            genes = pd.read_csv(genes, sep='\t', header=None)[0].tolist()

            #pandas load, transform to dict and remove the ones associated to non-coding
            regions = pd.read_csv(regions, sep='\t', header=None)
            regions = regions[regions[6] == 'protein-coding']

            if ARGS['regtype']:
                regions = regions[regions[4] == ARGS['regtype']]

            regions = regions.groupby(3)[5].apply(set).to_dict()

            print('genes =', len(genes), 'regions =', len(regions))
            print('***************************************************************')


            ng, kg = hypergeom_genome(genes, regions, ARGS['totalCodingGenes'])


            print('----------------------------------------------------------------')
            print('Hypergeometric test over genes for enrichment vs bckgd')

            hyper_background(genes, allpeaks, ng, kg)
            print('----------------------------------------------------------------')
            print('\n')