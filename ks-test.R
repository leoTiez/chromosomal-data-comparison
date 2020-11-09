###############################################################################
# Load/install libraries and files
###############################################################################
suppressMessages(library("optparse"))

option_list = list(
  make_option(c('-o', '--original_file'), type='character', default=NULL, 
    help='Path to original file', metavar='character'),
  make_option(c('-c', '--comparison_file'), type='character', default=NULL, 
              help='Path to comparison file'),
  make_option('--bin_size', type='numeric', default=0.5, 
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
  install.packages('ggplot2')
}

suppressMessages(library('rtracklayer'))
suppressMessages(library('devtools'))
suppressMessages(library('ggplot2'))

org_df <- import.bw(opt$original_file)
ref_df <- import.bw(opt$comparison_file)

org_df <- import.bw(org_file)
ref_df <- import.bw(comp_file)

org_scores <- score(org_df)
ref_scores <- score(ref_df)

org_scores <- (org_scores - mean(org_scores)) / sd(org_scores)
ref_scores <- (ref_scores - mean(ref_scores)) / sd(ref_scores)

print('############### KS test w/out binning')
ks.test(org_scores, ref_scores)

# KS based on previously calculated CDF
uv4_cdf <- ecdf(org_scores)(seq(min(org_scores), max(org_scores), by=opt$bin_size))
uv3b_cdf <- ecdf(ref_scores)(seq(min(ref_scores), max(ref_scores), by=opt$bin_size))

lower <- max(c(min(org_scores), min(ref_scores)))
upper <- min(c(max(org_scores), max(ref_scores)))

plot(
  seq(lower, upper, by=0.01),
  ecdf(org_scores)(seq(lower, upper, by=0.01)),
  type='l',
  xlab='ChIP Value',
  ylab='CDF',
  col='green'
)
lines(
  seq(lower, upper, by=0.01),
  ecdf(ref_scores)(seq(lower, upper, by=0.01)),
  col='blue'
)

plot(density(org_scores), col='green')
lines(density(ref_scores), col='blue')

print('############### KS test w/ binning')
ks.test(uv4_cdf, uv3b_cdf)
