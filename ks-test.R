###############################################################################
# Normalisation functions
###############################################################################

centre <- function(data){
  return((data - mean(data)) / sd(data))
}

remap <- function(data){
  n_data <- data - min(data)
  return(n_data / max(n_data))
}

###############################################################################
# Command line argument parser
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
              help='Install libraries if passed'),
  make_option(c('-n', '--normalise'), type='character', default='remap', 
              help='The normalisation method that should be used.
              Possible options are: remap (default) and center.
              remap remaps the all values between 0 and 1 with 0 the 
              minimum and 1 the maximum 
              (use this when using the KS-test with binning);
              center centres the data to zero mean and a std of 1',
              metavar='character')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

###############################################################################
# Load/install libraries and files
###############################################################################

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

###############################################################################
# Main
###############################################################################

org_df <- import.bw(opt$original_file)
ref_df <- import.bw(opt$comparison_file)

org_scores <- score(org_df)
ref_scores <- score(ref_df)

if(tolower(opt$normalise) == 'remap'){
  print('############### Use remap normalisation')
  org_scores <- remap(org_scores)
  ref_scores <- remap(ref_scores)
} else if(tolower(opt$normalise) == 'center'){
   print('############### Use center normalisation')
   org_scores <- centre(org_scores)
   ref_scores <- centre(ref_scores)
} else {
  stop('Invalid normalisation method. Pass remap or center.')
}

print('############### KS test w/out binning')
ks.test(org_scores, ref_scores)

# KS based on previously calculated CDF
uv4_cdf <- ecdf(org_scores)(seq(min(org_scores), max(org_scores), by=opt$bin_size))
uv3b_cdf <- ecdf(ref_scores)(seq(min(ref_scores), max(ref_scores), by=opt$bin_size))

plot(
  seq(min(uv4_cdf), max(uv4_cdf), by=opt$bin_size),
  ecdf(org_scores)(seq(min(uv4_cdf), max(uv4_cdf), by=opt$bin_size)),
  type='l',
  xlab='ChIP Value',
  ylab='CDF',
  col='green'
)
lines(
  seq(min(uv3b_cdf), max(uv3b_cdf), by=opt$bin_size),
  ecdf(ref_scores)(seq(min(uv3b_cdf), max(uv3b_cdf), by=opt$bin_size)),
  col='blue'
)

plot(density(org_scores), col='green')
lines(density(ref_scores), col='blue')

print('############### KS test w/ binning')
sprintf('############### Number of bins orginial data: %i', length(uv4_cdf))
sprintf('############### Number of bins reference data: %i', length(uv3b_cdf))
ks.test(uv4_cdf, uv3b_cdf)
