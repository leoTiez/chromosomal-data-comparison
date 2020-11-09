###############################################################################
# Load/install libraries and files
###############################################################################
suppressMessages(library("optparse"))

option_list = list(
  make_option(c('-o', '--original_file'), type='character', default=NULL, 
    help='Path to original file', metavar='character'),
  make_option(c('-c', '--comparison_file'), type='character', default=NULL, 
              help='Path to comparison file'),
  make_option('--bin_size', type='numeric', default=0.1, 
              help='Bin size'),
  make_option('--install', action='store_true', default=FALSE, 
              help='Install libraries if passed')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(opt$install == TRUE){
  source('http://www.bioconductor.org/biocLite.R')
  if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
  
  BiocManager::install('rtracklayer')
  install.packages('optparse')
}

suppressMessages(library('rtracklayer'))
suppressMessages(library('devtools'))

org_path <- opt$original_file
ref_path <- opt$comparison_file

org_df <- import.bw(org_path)
ref_df <- import.bw(ref_path)
org_scores <- (score(org_df) - mean(score(org_df))) / var(score(org_df))
ref_scores <- (score(ref_df) - mean(score(ref_df))) / var(score(ref_df))
print('############### KS test w/out binning')
ks.test(org_scores, ref_scores)

# KS based on previously calculated CDF
uv4_cdf <- ecdf(org_scores)(seq(min(org_scores), max(org_scores), by=opt$bin_size))
uv3b_cdf <- ecdf(ref_scores)(seq(min(ref_scores), max(ref_scores), by=opt$bin_size))
print('############### KS test w/ binning')
ks.test(uv4_cdf, uv3b_cdf)
