# Bait-ER

A fully Bayesian approach to estimate selection coefficients from Evolve-and-Resequence time series data.


## Version 

0.1 (December 2019)


## Citation

Soonish

## Download and compile Bait-ER

First of all, we need to download all the necessary files to compile and run Bait-ER. These files are in this GitHub repository. You can download these files manually or by using Git on the Ubuntu terminal:


```
git clone https://github.com/mrborges23/Bait-ER.git Bait-ER
```

Before compiling Bait-ER, we need to make sure that some packages are installed in your machine. Bait-ER uses `armadillo`, which requires packages like `LAPACK`, `boost` and `BLAS`. You can simply run the following commands (same applies to Mac OS X but using `brew`):


```
sudo apt-get install liblapack-dev
sudo apt-get install libblas-dev
sudo apt-get install libboost-dev
sudo apt-get install libarmadillo-dev
```

Bait-ER is coded in C++ and needs to be compiled. To do that, you can type the following command in your terminal:


```
g++ baiter.cpp -o baiter -O2 -std=c++14 -llapack -lblas -larmadillo
```

You may get the message that some libraries (e.g. `libhdf5.so`) are missing. You can solve this problem by finding the location of the missing files and adding them to the `$LD_LIBRARY_PATH path`. An example:


```
locate libhdf5.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path/to/missing/file
```

If the compilation proceeds without errors, you should have created an executable called `baiter` in your working directory. Before running Bait-ER, let us first explain its input files: the sync and control files.


## Sync file

Sync files are specified in Kofler et al. (2011). Sync files contain <img src="https://render.githubusercontent.com/render/math?math=3+n">  columns with column 1 indicating the chromosome (reference contig), column 2 the indicating position (in the reference contig) and column 3 indicating the reference allele. The following <img src="https://render.githubusercontent.com/render/math?math=n">  columns have the sync entries for each replicate and time point in the form of `A:T:C:G:N:deletion` counts. Sync files originally do not have a header, but headers are accepted when specified in the control file (value `1`).

An example sync file with 99 loci is provided to test Bait-ER. This data was taken from Burke et al. (2014):

```
chr  pos  ref  T0_R1           T0_R2           T0_R3           T0_R4	
2L   59   T    0:64:12:0:0:0   0:68:20:0:0:0   0:82:32:0:0:0   0:64:20:0:0:0
2L   62   A    70:0:13:0:0:0   69:0:21:0:0:0   90:0:32:0:0:0   67:0:22:0:0:0
2L   153  A    112:0:28:0:0:0  108:0:26:0:0:0  153:0:37:0:0:0  79:0:28:0:0:0
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
|`Prior_parameters`|The prior parameters <img src="https://render.githubusercontent.com/render/math?math=\alpha"> and <img src="https://render.githubusercontent.com/render/math?math=\beta"> of a gamma distribution. These are the parameters of the gamma distribution that models the distribution of <img src="https://render.githubusercontent.com/render/math?math=\sigma">. | These parameters should be small enough. Simulations conducted by us showed that values around 0.001 are unlike to overinfluence the posterior for the generality of E&R experimental designs. |
|`Output_file`|The name of the output file. ||


## Running Bait-ER

Once you have your sync and control files on the same folder as the `baiter` executable, you can open the terminal and run the executable followed by the name of the control file:

```
./baiter baiter.cf
```

Bait-ER will immediately output some information. Confirm that this information conforms to your control and sync file. To make sure Bait-ER is running, you should get the message `Bait-ER has started!`. When the analyses are done, you should get the message `Bait-ER has finished!`. Notice, that Bait-ER periodically writes of the `output_file`.


## Output file

The output file has information about the posterior of <img src="https://render.githubusercontent.com/render/math?math=\sigma"> per locus. The output file has six columns with the following order: chromosome, position, reference, average <img src="https://render.githubusercontent.com/render/math?math=\sigma">, log Bayes factor (for the hypothesis that <img src="https://render.githubusercontent.com/render/math?math=\sigma"> is different from 0), and the posterior values of <img src="https://render.githubusercontent.com/render/math?math=\alpha"> and <img src="https://render.githubusercontent.com/render/math?math=\beta"> (also known as shape and rate parameters). The reference column reports the reference allele for which <img src="https://render.githubusercontent.com/render/math?math=\sigma"> was calculated, not the reference allele indicated on the `sync_file`.

The output for the first three loci of the example.txt data is the following:

```
chromosome position  reference  sigma      logBF      alpha         beta
2L         59        T:C        -0.000202  -0.044448  15511.266710  15514.402451
2L         62        A:C        -0.000506  -0.105363  15663.110987  15671.040724
2L         153       A:C        -0.003621  -0.757280  16573.253165  16633.477700
```

<img src="https://render.githubusercontent.com/render/math?math=\alpha"> and <img src="https://render.githubusercontent.com/render/math?math=\beta"> can be used to estimate other quantities of interest regarding the posterior distribution of <img src="https://render.githubusercontent.com/render/math?math=\sigma"> (or more correctly, <img src="https://render.githubusercontent.com/render/math?math=1+\sigma">, the fitness distribution): quantiles, credibility intervals, etc. To do so, one can use the `qgamma` function in `R`. For example, if we want a 95% credibility interval for <img src="https://render.githubusercontent.com/render/math?math=\sigma"> at position 59 of chromosome 2L, one can simply type in `R`:

```R
shape <- 15511.266710                  # shape parameter 
rate  <- 15514.402451                  # rate parameter
qgamma(0.05,shape=shape,rate=rate)-1   # lower bound
[1] -0.01336969
qgamma(0.95,shape=shape,rate=rate)-1   # upper bound
[1] 0.01303874
```

<img src="https://render.githubusercontent.com/render/math?math=\sigma"> varies between [-0.013,0.013] with 95% probability, which includes 0 (i.e. neutral evolution). This result is in line with the logBF on the output table, which by being close to zero, is not suggesting that this loci constitutes a target of selection. The absolute value of the logBFs can be used to conclude whether a single locus in evolves under neutrality or selection, just like the log(p-value) used to build the standard Manhattan plots.

Since the E&R studies include thousands to millions of loci, we need to be a little bit more stringent about the BFs thresholds to use to select targets of selection. Plese read next section for more information on how to properly correct BFs with Bait-ER.

## A note on BFs thresholds

To infer the genomic response to adaptation, Bait-ER employs multiple tests for a set of thousands of millions of loci throughout the genome. Multiple testing is prone to false positives, and we advise the users of Bait-ER to correct their BFs with proper multiple testing correction strategies. In particular, we recommend the strategy proposed by Wellcome Trust Case Control Consortium (2007). Their correction requires that we have an idea of the fraction <img src="https://render.githubusercontent.com/render/math?math=F"> of loci under selection. Then, we set an aimed power of the test <img src="https://render.githubusercontent.com/render/math?math=P"> and the odds <img src="https://render.githubusercontent.com/render/math?math=O"> in favor of a true positive to a false positive. These three quantities can be translated into a BF that serves as a threshold for picking targets of selection:

<img src="https://render.githubusercontent.com/render/math?math=BF_T=\log \left( \frac{O}{PF}-1 \right)">

If we consider that an E&R experiment should have 0.1% expected loci under selection, and set an aimed power of 0.5 and odds favoring true positives over false positives of 10, we calculate a corrected BF of 9.9. This threshold is considerably more restrictive than the typical thresholds used in standard Bayesian applications, where BFs higher than 3.4 represent already very strong evidence in favor of the tested hypothesis. 


## References

* Burke, Liti and Long (2014) Standing genetic variation drives repeatable experimental evolution in outcrossing populations of saccharomyces cerevisiae. Molecular Biology and Evolution 31(12): 3228–3239

* Jónás, Taus, Kosiol, Schlötterer and Futschik (2016) Estimating the effective population size from temporal allele frequency changes in experimental evolution. Genetics 204(2):723-735

* Kofler, Pandey, Schlötterer (2011) PoPoolation2: identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq). Bioinformatics 27(24):3435–3436

* Wellcome Trust Case Control Consortium (2007) Genome-wide association study of 14,000 cases of seven common diseases and 3,000 shared controls. Nature 447(7145):661-78


## Questions and bug reporting

Please use **Issues** to report possible bugs, suggest enhancement features, or if you need help using Bait-ER. If you have more theoretical or biological questions, you can directly contact Rui Borges (ruiborges23@gmail.com) or Carolina Barata (cdcbrb@st-andrews.ac.uk).


## License

This program is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software. See the GNU General Public License (http://www.gnu.org/licenses/) for more details.
