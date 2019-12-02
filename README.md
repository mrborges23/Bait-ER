# Bait-ER

A fully Bayesian approach to estimate selection coefficients from Evolve-and-Resequence allele trajectories.


## Download and compile Bait-ER

First of all, we need to download all the necessary files to compile and run Bait-ER. These files are in this GitHub repository. You can download these files manually or by using Git on the Ubuntu terminal:


```
git clone https://github.com/mrborges23/Bait-ER.git Bait-ER
```

Before compiling Bait-ER, we need to make sure that some packages are installed in your machine. Bait-ER uses `armadillo`, which requires packages like `LAPACK`, `boost` and `BLAS`. You can simply run the following commands:


```
sudo apt-get install liblapack-dev
sudo apt-get install libblas-dev
sudo apt-get install libboost-dev
sudo apt-get install libarmadillo-dev
```

Bait-ER is coded in C++ and needs to be compiled. To do that, you can type the following command in your terminal:


```
g++ main.cpp -o baiter -O2 -std=c++14 -llapack -lblas -larmadillo
```

You may get the error that some libraries (e.g. `libhdf5.so`) are missing. You can solve this problem by finding the location of the missing files and adding them to the `$LD_LIBRARY_PATH path`. An example:


```
locate libhdf5.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path/to/missing/file
```

If the compilation proceeds without errors, you should have created an executable called `baiter` in your working directory. Before running Bait-ER, let us first explain its input files: the sync and control files.


## Sync file

Sync files are specified in Kofler et al. (2011). Sync files contain 3 + n columns with column 1 indicating the chromosome (reference contig), column 2 the indicating position (in the reference contig) and column 3 indicating the reference allele. The following n columns have the sync entries for each replicate and time point in the form of `A:T:C:G:N:deletion` counts. Sync files originally do not have a header but headers are accepted when specified in the control file (value `1`).

```
chr  pos  ref R1_T0          R2_T0          R3_T0          R4_T0
chrI 3049 A   233:70:0:0:0:0 233:70:0:0:0:0 233:70:0:0:0:0 33:70:0:0:0:0
chrI 3056 C   67:0:0:236:0:0 67:0:0:236:0:0 67:0:0:236:0:0 67:0:0:236:0:0
chrI 3141 C   0:84:0:282:0:0 0:84:0:282:0:0 0:84:0:282:0:0 0:84:0:282:0:0
```

## Control file

The control file includes all the necessary parameter to run Bait-ER:

|Parameter|Description|Notes|
|---      |---        |---  |
|`Sync_file`|The name of the sync file that we want to analyze.||
|`Header`|This parameter indicates whether the sync file has (value `1`) or has not (value `0`) a header.||
|`Columns_order`| The way columns are organized in the sync file. Columns can be organized in time point vs. replicate (value `0`) or replicate vs. time point (value `1`). | Example. Consider that we have an experimental design with two-time points and two replicates.  In order `0`, Bait-ER expects the columns to be organized like this: T1R1, T1R2, T2R1 and T2R2. In order `1`, Bait-ER expects the columns to be organized like this: T1R1, T2R1, T1R2 and T2R2. |
|`Number_replicates`|The number of replicates| |
|`Time_points` |The vector of time points. Their values should be separated by commas. | Be careful for not introducing any space or other characters after or before the comma.|
|`Number_loci` |The number of loci in the sync file.| |
|`Population_size` |The effective population size. |This can be estimated by other methods (e.g. Jónás et al. (2016)). |
|`Prior_parameters`|The prior parameters alpha and beta of a gamma distribution. These are the parameters of the gamma distribution that models the distribution sigma. | These parameters should be small enough. Simulations conducted by us showed that values smaller or equal than 0.001 are unlike to over influence the posterior for the generality of E&R experimental designs. |
|`Output_file`|The name of the output file. ||


## Running Bait-ER

Once you have your sync and control files on the same folder as the baiter executable, you can open the terminal and run the executable:

```
./baiter 
```

Bait-ER will immediately output some information. Confirm that this information conforms your control and sync file. To make sure Bait-ER is running, you should get the message `Bait-ER has started!`. When the analyses are done, you should get the message `Bait-ER has finished!`. Notice, that Bait-ER periodically writes of the `output_file`. To check whether your analysis is frozen, you can simply check if the output has been updating. 


## Output file

The output file has information about the posterior of sigma per locus. The output file has six columns with the following order: chromosome, position, average sigma, log Bayes factor (for the hypothesis that sigma is different from 0), and the posterior values of alpha and beta (also know as shape and rate parameters). 

The output for the first three loci of the example.txt data is the following:

```
chromosome position sigma      logBFs   alpha   beta
chrI       3049     -0.0129154 -1.57228 5229.56 5297.98
chrI       3056     -0.0184936 -2.37067 5312.73 5412.84
chrI       3141     -0.0092678 -1.08498 5035.45 5082.55
```

Alpha and beta can be used to estimate other quantities of interest regarding the posterior distribution of sigma (or more correctly, 1 + sigma, the distribution of fitness): quantiles, intervals of confidence, etc. To do so, one can use the gamma function in `R`. For example, if we want a 95% credibility interval for sigma at position 3049, one can simply type in R:

```
shape <- 5229.56 # shape parameter 
rate  <- 5297.98 # rate parameter
qgamma(0.05,shape=shape,rate=rate)-1  # lower bound
[1] -0.03525825
qgamma(0.95,shape=shape,rate=rate)-1  # upper bound
[1] 0.009644153
```

One can conclude that sigma is between [-0.035,0.009] with 95% probability, which includes 0. This result is in line with the logBF on the table, that for being close to zero, is not suggesting that this position constitutes a target of selection. 

## References

* Jónás, Taus, Kosiol, Schlötterer and Futschik (2016) Estimating the effective population size from temporal allele frequency changes in experimental evolution. Genetics 204(2):723-735

* Kofler, Pandey, Schlötterer (2011) PoPoolation2: identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq). Bioinformatics 27(24):3435–3436


## Version 

0.0 (December 2019)


## Citation

Soonish


## Questions and bug reporting

Please use **Issues** to report possible bugs, suggest enhancement features or if you need help using Bait-ER. If you have more theoretical or biological questions, you can directly contact Rui Borges (ruiborges23@gmail.com) or Carolina Barata (cdcbrb@st-andrews.ac.uk).


## License

This program is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software. See the GNU General Public License (http://www.gnu.org/licenses/) for more details.
