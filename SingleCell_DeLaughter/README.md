
# Scripts for Single Cell data processing:


### Mouse reference genome 
was obtained from the following link: 
http://ftp.ensembl.org/pub/current_gtf/mus_musculus/



### experiments.xlsx 
includes updates on sample information



### library_builder.R 
builds the transcripts to genes library for next step



### t2gene.R 
merges transcripts into genes, generates TPMs and Est. counts for each gene



### high_var_gene_filter.R 
applies both linear regression and sigmoid fit to select highly varible genes 



### mds_plots.R 
does multidimensional scaling on normalized dataset, and plot with perfect genes markers for each cell type



### total_reads_counts.py 
counts total number of reads for each cell



### cv^2_vs_mean_plot_setup.txt
setting and parameters used for fitting



### only normal cells are included
