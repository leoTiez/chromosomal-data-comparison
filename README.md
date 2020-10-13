# Comparison of chromosomal data

This program compares your own data with the
a reference data set that is produced by
Wyrick et al. 2016

## Requirements 
You need to have the following requirements installed
- Python >= 3.6
- pip3

## Installation
To install the necessary python libraries, run
```bash
python3 -m pip install -r requirements
```

## Execution
Execute the file via
                                                        
```
python3 main.py --bed=file/path -i=file/path [-i=more/file/paths] [-n="YOUR NAME" [-name="More names matching input"]] [--smoothing=200 [--smoothing=None]] [--thresh=0.95] [--save_plot] [--save_prefix='Your_prefix'] [--num_lags=150]
```
- `--bed`: Required: Path to the bed file.
- `--input_data``-i`: Required. Path to input file. Add more inputs via more `-i=next/file` or `--input_data=next/file`
- `--name` `-n`: Optional. Names for creating better plots. Otherwise numbers are used. Add more names via `--name="Next Name"` or `-n="Next Name"`. Note that the number of names (if passed) must match the number of inputs
- `--smoothing`: Optional. Number of values that are used for smoothing. None per default (not applied). Add more smoothing factors via `--smoothing=200` where 200 is replaced by your value. Note that the number of smoothing factors (if passed) must match the number of input data.
- `--thresh`: Optional. Threshold that is set for how much percent of the data should be matched by distance measure to be defined as a similar signal. Default is 0.9.
- `--save_plot`: Optional. Flag is set when the plots are to be saved
- `--save_prefix`: Optional. Set a prefix for the saved plots as an identifier. Not used if `--save_plot` is not set
- `--num_lags`: Optional. Maximal number of values that the signal is shifted for the MSE. Default is 100.

All paths must be relative to your current directory.