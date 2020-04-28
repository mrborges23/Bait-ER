# Bait-ER

A fully Bayesian approach to estimate selection coefficients from Evolve-and-Resequence time series data.


## Version 

0.1 (December 2019)


## Citation

Soonish

## Download and compile Bait-ER

First of all, we need to download all the necessary files to compile and run Bait-ER; these files are in this GitHub repository. You can download these files manually or by using Git on the Ubuntu terminal:


```
git clone https://github.com/mrborges23/Bait-ER.git Bait-ER
```

Before compiling Bait-ER, we need to make sure that some packages are installed on your computer. Bait-ER uses `armadillo`, which requires the `LAPACK`, `boost` and `BLAS` packages. The following commands install these packages (same applies to Mac OS X, but using `brew` instead of `sudo apt-get`):


```
sudo apt-get install liblapack-dev
sudo apt-get install libblas-dev
sudo apt-get install libboost-dev
sudo apt-get install libarmadillo-dev
```

Bait-ER is coded in C++ and needs to be compiled; this step will produce an executable. We use the `g++` compiler, but any other can be used:

```
g++ baiter.cpp -o baiter -O2 -std=c++14 -llapack -lblas -larmadillo
```

You should have created an executable called `baiter` in your working directory. If you experience erros during the compilation step (e.g., libraries missing), please read the section **Compilation errors and how to solve them**. Before running Bait-ER, let us first explain its input files: the sync and control files.


## Sync file

Sync files are described in Kofler et al. (2011). Sync files contain 3+n columns: column 1 indicates the chromosome (reference contig); column 2 registers the position (in the reference contig); column 3 shows the reference allele. The following n columns include the observed counts for each replicate and time point of the experiment, following the format  `A:T:C:G:N:deletion`. Sync files originally do not have a header, but headers are accepted when specified in the control file (`Header` with value `1`).

An example sync file with 100 loci is provided to test Bait-ER; this data was taken from Barghi et al. (2019) and regard an E&R experiment in hot-adapted populations of *Drosophila simulans* (10 replicates and 7 times points sampled at each 10 generations):

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
|`Sync_file`|The name of the sync file that will be analyzed.||
|`Header`|This parameter indicates whether the sync file has (value `1`) or has not (value `0`) a header.||
|`Columns_order`| Indicates the way columns are organized in the sync file. Columns can be organized by time point vs. replicate (value `0`) or by replicate vs. time point (value `1`). | Example: Consider that we have an experimental design with two-time points and two replicates.  In order `0`, Bait-ER expects the columns organized like this: T1R1, T1R2, T2R1, and T2R2. In order `1`, Bait-ER expects the columns organized like this: T1R1, T2R1, T1R2, and T2R2. |
|`Number_replicates`|The number of replicates. | |
|`Time_points` | The vector of time points in generations: commas should separate these values. | Be careful for not introducing any space or other characters after or before the comma. |
|`Number_loci` | The number of loci in the sync file. | |
|`Population_size` |The effective population size. | The effective population size can be estimated by other methods: e.g., Jónás et al. (2016). |
|`Prior_parameters`| The vector of prior parameters: commas should separate these values. The prior parameters <img src="https://render.githubusercontent.com/render/math?math=\alpha"> and <img src="https://render.githubusercontent.com/render/math?math=\beta"> (also known as shape and rate parameters) set a prior gamma distribution on <img src="https://render.githubusercontent.com/render/math?math=\sigma">. | Simulations conducted by us showed that values around 0.001 are unlike to overinfluence the posterior, for the generality of E&R experimental designs. |
|`Output_file`| The name of the output file. ||


## Running Bait-ER

Place your sync and control files on the same folder that has the Bait-ER executable. Then, open the terminal and run the executable `baiter` followed by the name of the control file:

```
./baiter baiter.cf
```

Bait-ER will immediately output some information. Confirm that this information conforms to your control and sync file. To make sure Bait-ER is running, you should get the message `Bait-ER has started!`. When the analyses terminate, you should get the message `Bait-ER has finished!`. Bait-ER periodically writes of the `output_file`.


## Output file

The output file has information about the posterior of <img src="https://render.githubusercontent.com/render/math?math=\sigma"> per locus. It has six columns with the following order: chromosome, position, reference, average <img src="https://render.githubusercontent.com/render/math?math=\sigma">, log Bayes factor (for the hypothesis that <img src="https://render.githubusercontent.com/render/math?math=\sigma"> is different from 0), and the posterior values of <img src="https://render.githubusercontent.com/render/math?math=\alpha"> and <img src="https://render.githubusercontent.com/render/math?math=\beta"> (also known as shape and rate parameters). The reference column reports the reference allele for which <img src="https://render.githubusercontent.com/render/math?math=\sigma"> was calculated, this is the allele that had the highest observed counts (not necessarily the major or rising allele). The alleles following `:` are the remaining observed alleles at that position. Bait-ER handels triallelic and fourallelic sites.

The output for the first three loci of the example.txt data is the following:

```
chromosome position  reference  sigma      logBF      alpha         beta
2L         59        T:C        -0.000202  -0.044448  15511.266710  15514.402451
2L         62        A:C        -0.000506  -0.105363  15663.110987  15671.040724
2L         153       A:C        -0.003621  -0.757280  16573.253165  16633.477700
```

The absolute value of the logBFs (fifth columns) can be used to conclude whether a single locus in evolves under neutrality or selection, just like the log(p-value) used to build the standard Manhattan plots. An example of the chromosome 2L of hot-adapted D. simulans populations follows (Barghi et al. 2019): the targets of selection are highlighted in red alongside with their positions in the chromosome 2L.

![Manhattan_plot](https://github.com/mrborges23/Bait-ER/blob/master/Manhattan_plot.jpeg)

Statistical significance is also assessed via the Bfs; however, the standard BF thresholds a quite relaxed (e.g., log(99)). Since the E&R studies include thousands to millions of loci, we need to be a little bit more stringent about the BFs thresholds to use to select targets of selection. Please read the section * A note on BFs thresholds* for more information on how to correct BFs with Bait-ER properly.

## Estimating other statistics of interest for the selection coefficients

Bait-ER exports the two statistics that are needed to perform inferences: the average of <img src="https://render.githubusercontent.com/render/math?math=\sigma">  and the log BFs (fourth and fifth columns). However, we also output the shape and rate (<img src="https://render.githubusercontent.com/render/math?math=\alpha"> and <img src="https://render.githubusercontent.com/render/math?math=\beta">) parameters of the posterior gamma distribution of <img src="https://render.githubusercontent.com/render/math?math=\sigma">, which can be used to calculate other quantities of interest (quantiles, credibility intervals, etc). The gamma distribution is defined over the fitness domain: i.e., 1 + <img src="https://render.githubusercontent.com/render/math?math=\sigma">, so be carefull to substract 1 whenever you want to report a statistic for <img src="https://render.githubusercontent.com/render/math?math=\sigma"> . 

To calculate additional statistics one can use the `qgamma` function in `R`. For example, if we want a 95% credible interval for <img src="https://render.githubusercontent.com/render/math?math=\sigma"> at position 59 of chromosome 2L (check the output file for the values of <img src="https://render.githubusercontent.com/render/math?math=\alpha"> and <img src="https://render.githubusercontent.com/render/math?math=\beta">  in this position), one can simply use these `R` commands:

```R
shape <- 15511.266710                  # shape parameter 
rate  <- 15514.402451                  # rate parameter
qgamma(0.05,shape=shape,rate=rate)-1   # lower bound
[1] -0.01336969
qgamma(0.95,shape=shape,rate=rate)-1   # upper bound
[1] 0.01303874
```

<img src="https://render.githubusercontent.com/render/math?math=\sigma"> varies between [-0.013,0.013] with 95% probability, which includes 0 (i.e. neutral evolution). This result is in line with the logBF on the output table, which by being close to zero, is not suggesting that this loci constitutes a target of selection. 

## A note on BFs thresholds

To infer the genomic response to adaptation, Bait-ER employs multiple tests for a set of thousands of millions of loci throughout the genome. Multiple testing is prone to false positives, and we advise the users of Bait-ER to correct their BFs with proper multiple testing correction strategies. In particular, we recommend the strategy proposed by Wellcome Trust Case Control Consortium (2007). Their correction requires that we have an idea of the fraction <img src="https://render.githubusercontent.com/render/math?math=F"> of loci under selection. Then, we set an aimed power of the test <img src="https://render.githubusercontent.com/render/math?math=P"> and the odds <img src="https://render.githubusercontent.com/render/math?math=O"> in favor of a true positive to a false positive. These three quantities can be translated into a BF that serves as a threshold for picking targets of selection:

<img src="https://render.githubusercontent.com/render/math?math=BF_T=\log \left( \frac{O}{PF}-1 \right)">

If we consider that an E&R experiment should have 0.1% expected loci under selection, and set an aimed power of 0.5 and odds favoring true positives over false positives of 10, we calculate a corrected BF of 9.9. This threshold is considerably more restrictive than the typical thresholds used in standard Bayesian applications, where BFs higher than 3.4 represent already very strong evidence in favor of the tested hypothesis. 

## Compilation errors and how to solve them

While compiling Bait-ER, you may get the message that some libraries (e.g. `libhdf5.so`) are missing. You can solve this issue by finding the location of the missing files and adding it to the `$LD_LIBRARY_PATH path`. An example:

```
locate libhdf5.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path/to/missing/file
```

## References

* Barghi, Tobler, Nolte, Jakšić, Mallard, Otte, Dolezal, Taus, Kofler, Schlötterer (2019) Genetic redundancy fuels polygenic adaptation in Drosophila. PLoS Biology 17(2):e3000128

* Jónás, Taus, Kosiol, Schlötterer and Futschik (2016) Estimating the effective population size from temporal allele frequency changes in experimental evolution. Genetics 204(2):723-735

* Kofler, Pandey, Schlötterer (2011) PoPoolation2: identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq). Bioinformatics 27(24):3435–3436

* Wellcome Trust Case Control Consortium (2007) Genome-wide association study of 14,000 cases of seven common diseases and 3,000 shared controls. Nature 447(7145):661-78


## Questions and bug reporting

Please use **Issues** to report possible bugs, suggest enhancement features, or if you need help using Bait-ER. If you have more theoretical or biological questions, you can directly contact Rui Borges (ruiborges23@gmail.com) or Carolina Barata (cdcbrb@st-andrews.ac.uk).


## License

This program is free software. You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software. See the GNU General Public License (http://www.gnu.org/licenses/) for more details.

