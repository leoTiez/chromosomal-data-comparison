# Comparison of chromosomal data

This program is designed to provide different metrics
to compare similarity or dissimilarity of ChIP-seq data. 
Applied measurements are thresholding (how large does a tolerated
error margin need to be in order to account for k% of the data); mean squared
error (MSE); and KS-testing. Note that if you want to compare the
ChIP-seq signal of the whole genome, the calculated p-values will
be (close to) zero. In this case, binning is recommended.


There are two scripts: one Python script, one Rscript. The Python script
is recommended to use as it provides more flexibility. However, the Rscript
gives the possibility to compare the outcome of the KS-test with 
another implementation to rule out numerical issues.

## Requirements 
You need to have the following requirements installed
- Python >= 3.6
- pip3

If you want to use the R scripts you need to install
- R >= 3.0.0

If you install R for the first time, run afterwards
```bash
sudo apt -y install libcurl4-gnutls-dev libxml2-dev libssl-dev
```
and install `devtools` with root privileges

```bash
sudo R
> packages.install('devtools')
```

## Installation
To install the necessary python libraries, run
```bash
python3 -m pip install -r requirements
```

## Python Execution
Execute the file via
                                                        
```
python3 main.py -i=file/path [-i=more/file/paths] [-n="YOUR NAME" [-name="More names matching input"]] [-bed="path/to/bed"] [--smoothing=200 [--smoothing=None]] [--thresh=0.95] [--save_plot] [--save_prefix='Your_prefix'] [--num_lags=150] [--norm=remap] [--num_bins=10]
```
- `--input_data``-i`: Required. Absolute path to input file. Add more inputs via more `-i=next/file` or `--input_data=next/file`
- `--name` `-n`: Optional. Names for creating better plots. Otherwise numbers are used. Add more names via `--name="Next Name"` or `-n="Next Name"`. Note that the number of names (if passed) must match the number of inputs
- `--bed`: Optional. Path to bed file. If passed, mean and std are computed for the segments defined in the file.
- `--smoothing`: Optional. Number of values that are used for smoothing. None per default (not applied). Add more smoothing factors via `--smoothing=200` where 200 is replaced by your value. Note that the number of smoothing factors (if passed) must match the number of input data.
- `--thresh`: Optional. Threshold that is set for how much percent of the data should be matched by distance measure to be defined as a similar signal. Default is 0.9.
- `--save_plot`: Optional. Flag is set when the plots are to be saved
- `--save_prefix`: Optional. Set a prefix for the saved plots as an identifier. Not used if `--save_plot` is not set
- `--num_lags`: Optional. Maximal number of values that the signal is shifted for the MSE. Default is 100.
- `--norm`: Optional. The normalisation method that is used. If none is passed, no normalisation is applied. Otherwise 
set it to `center` for zero mean and std of one, or `remap` for rescaling values between 0 and 1
- `--num_bins`: the number of bins that are used for the KS-test. If none is passed, the `doane` binning method is applied
(see numpy documentation). Otherwise you can pass any positive integer value.

All paths must be relative to your current directory.

## R execution
The R script creates two plots, the CDF and the probability density function
and saves them in a PDF in your working directory. The R script does not accept a bed
file at the moment to aggregate mean and std over several segments.
If you want to run the R script to compare the
outcome of the KS-test, execute the following command

```bash
Rscript -o [--original_file] file/to/bwFile1 -c [--comparison-file] file/to/bwFile2 [--bin_size 0.1] [--install] [-n [--normalise] remap]
```
 where 
 - `-o` or `--original_file`: Required. Determines the path to your original bigwig file
 - `-c` or `--comparison_file` Required. Determines the path to the bigwig file you want to compare your original file with
 - `--bin_size`: Optional. Sets the range per bin
 - `--install`: Optional. If set, all necessary libraries are installed
 - `-n` or `--normalise`: Optional. If not passed, no normalisation method is applied. Otherwise, choose between `remap` 
 (rescaling the data between 0 and 1), or `center` (zero mean and std of 1)