#!/usr/bin/env python

#AUTHOR: Nicole Ferraro
#BMI 217 Project, Winter 2016-2017
#PURPOSE: Concatenate various data types from downloaded files from the CCLE site into one input file for use as training in downstream models
#No normalization is done in this script, done with model training

import csv
import numpy as np
import pandas as pd
import time

#Assumes data exists one folder up in directory called CCLE, change as needed
exp_data = pd.read_csv('../CCLE/CCLE_Expression_Entrez_2012-10-18.res', sep='\t', index_col=False)
        
cnv_data = pd.read_csv('../CCLE/CCLE_copynumber_byGene_2013-12-03.txt', sep='\t')
mut_data = pd.read_csv('../CCLE/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf', sep='\t', index_col=False)
drug_data = pd.read_csv('../CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv', sep=',')

#Only extract genes with both expression and CNV data, and shared cell lines
genes1 = cnv_data['SYMBOL'].unique().tolist()
genes2 = exp_data['Description'].unique().tolist()
overlap_genes = set(genes1).intersection(genes2)
mut_genes = mut_data['Hugo_Symbol'].unique()
cells = drug_data['CCLE Cell Line Name'].unique()

new_cols = cells.tolist().insert(0,'Gene')

#Sort each file by cell line (treat as samples)
#Get overlap of cell lines
#Combine data where row is gene, repeated 2-6 times, to get above three types of data (three types of mutation data)
#Across all cell lines (columns)
#Long run time (~ 1 hr)

start = True
for gene in overlap_genes:
    #These are the six possible lines for each gene
    #All genes will have expr and cnv
    expr = [gene, 'Expression']
    cnv = [gene, 'CNV']
    mut = [gene, 'MutLOF']
    mut_nn = [gene, 'nnMS']
    mut_utr = [gene, 'UTR']
    gene_mut = False
    if start == True:
      inc_cells = ['Gene', 'Type']
    #Loop over each cell line
    for cell in cells:
        if cell not in exp_data.columns or cell not in cnv_data.columns:
            continue
        if start == True:
          inc_cells.append(cell)
        #Get expression and cnv data
        ind1 = np.where(exp_data['Description'] == gene)
        gene_exp = exp_data[cell].iloc[ind1]
        ind2 = np.where(cnv_data['SYMBOL'] == gene)
        gene_cnv = cnv_data[cell].iloc[ind2]
        #Check if gene is one of the ones with sequencing and mutation data
        if gene in mut_genes:
          inds3 = np.where(mut_data['Hugo_Symbol'] == gene)
          inds4 = np.where(mut_data['Tumor_Sample_Barcode'] == cell)
          inds = np.intersect1d(inds3[0], inds4[0])
          if len(inds) == 0:
             gene_mut = False
          else:
              gene_mut = True
              mutLOF = 0
              nnMS = 0
              utr = 0
              #Check if one of the mutation types of interest is present in each sample for current gene
              for ind in inds:
                 mut_type = mut_data['Variant_Classification'].iloc[ind]
                 if 'Nonsense' in mut_type or 'Frame_Shift' in mut_type or 'Splice_Site' in mut_type:
                   mutLOF = 1
                 elif 'Missense' in mut_type or 'Intron' in mut_type or 'In_Frame' in mut_type:
                   nnMS = 1
                 elif 'UTR' in cell:
                   utr = 1
        if len(gene_exp) > 0 and len(gene_cnv) > 0:
          expr.append(gene_exp.iloc[0])
          cnv.append(gene_cnv.iloc[0])
          if gene_mut == True:
            mut.append(mutLOF)
            mut_nn.append(nnMS)
            mut_utr.append(utr)
    #Write out new data in current directory row by row
    with open('CCLE_data_format_mutation.txt', 'a') as f:
      writer = csv.writer(f, delimiter='\t')
      if start == True:
        writer.writerow(inc_cells)
        start = False
      writer.writerow(expr)
      writer.writerow(cnv)
      if gene_mut == True:
        writer.writerow(mut)
        writer.writerow(mut_nn)
        writer.writerow(mut_utr)

