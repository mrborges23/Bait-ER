# Bait-ER

A fully Bayesian approach to estimate selection coefficietns from Evolve-and-Resequence allele trajectories.

## Download and compile Bait-ER

First of all, we need to download all the necessary files to compile and run Bait-ER. These file are in this GitHub repository. You can download these files manually or by using Git on the Ubuntu terminal:


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

Bait-ER is coded in C++ and needs to be compiled. To do that, you can type the following comman in your terminal:


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



## Control file

The control file includes all the necessary parameter to run Bait-ER:

|Parameter|Description|Notes|
|---      |---        |---  |
|`Sync_file`|The name of the sync file that we want to analyze.||
|`Header`|This parameter indicates whether the sync file has (value 1) or has not (value 0) a header.||
|`Columns_order`| The way columns are organized in the sync file. Columns can be organized in time point vs. replicate (value 0) or replicate vs. time point (value 1). |
Example. Consider that we have an experimental design with two time points and two replicates.  In order 0, we expect the columns to dispose like this: T1R1, T1R2, T2R1 and T2R2. In order 1 we expect the columns to be dispose like this: T1R1, T2R1, T1R2 and T2R2. |
|`N_replicates`|The number of replicates| |
|`Time_points` |The vector of time points. Their values should be separated by commas. | Be careful for not introducing any space or other character after or before the comma.|
|`N_loci` |The number of loci in the sync file.| |
|`Population_size` |The effective population size. |This can be estimated by other methods: e.g. Agnes paper here. |
|`Prior_parameters`|The prior parameters alpha and beta of a gamma distribution. These are the parameters of the gamma distribution that models the distribution sigma. | These parameters should be small enough. Simulations conducted by us showed that values smaller or equal than 0.001 are unlike to overdominate the posterior for the generality of E&R experimental designs. |
|`Output_file`|The name of the output file. ||






## Running Bait-ER

Once you have your sync and control file, onling with the baiter executable, in the same folder you can now open the termal and run the baiter executable:

```
./baiter 
```

Bait-ER will immediatly output some information. Confirm that this information conforms your control file. If not, check whether your control file is correctly filled. 

A will be printer periodically, so you know how long Bait-ER will take to finish the analyses. However, the output file is updated constantly. If for some reason, your runs are stopped, you can always start by the last site that was printed to the update file.  




This codes can by used by following the main.R script. You can used Bait-ER to estimate selection coefficients from both simulated allele trajectories or a sync file. Toy data is provided!

Please contact Rui Borges (ruiborges23@gmail.com) or Carolina Barata (cdcbrb@st-andrews.ac.uk) if you have any questions regarding the model/code/output.





## Version 
0.0 (December 2019)


## Citation

Soonish


## Questions and bug reporting

Please use Issues to report possible bugs, enhancement features, help or simply questions reagarding Bait-ER. If zour question have nothing to do zou with Bait-ER but more with the model, data and output, zou can contact us directly: Rui Borges (ruiborges23@gmail.com) or Carolina Barata (cdcbrb@st-andrews.ac.uk).


## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software. See the GNU General Public License (http://www.gnu.org/licenses/) for more details.
