## New Order of Events

* Need to install Qiime 
	* use Pat's Installation
	* Install locally without sudo
	* Ask administrators of axiom to install
	
**For this I have decided to look at using Anaconda and then using this to install qiime**

So what did I do?

```bash
# Download the respective Anaconda Python 2.7 package for Linux
# Put into home directory in axiom

bash Anaconda2-2.5.0-Linux-x86_64.sh

```
From here I created a few alias' for a number of commands in my .bash_profile.
The first was to alias python2.7, second I aliased conda
Second I added anaconda2 to my PATH just to be on the safe side.

Next I followed instructions laid out by Jorge-C:

```bash
conda create -n qiime190rc2 python pip numpy matplotlib IPython future natsort scipy pandas scikit-bio gdata
```
I had to move into the qiime190rc2 bin directory to run the next command for me it was:
/home/marcsze/bin/anaconda2/envs/qiime190rc2/bin

```bash
source activate qiime190rc2
pip install https://github.com/biocore/qiime/archive/1.9.0-rc2.tar.gz
```
To make sure that this installed okay I ran a test as recommended by the Qiime website

```bash
print_qiime_config.py -t
```
The result returned *Ran 9 tests in 0.029s OK*.  This works regardless of the directory I am in so it seems to have worked.

The version that I have is missing a few dependencies (not sure if I need them though)
The following are what is missing: h5py, sortmerna, sumaclust

## Things to do...

* Decided not to use Qiime
	* Done

* Need to reorganize repository so that everything runs automatically
	* Done
	
* Need to update processing so all 454 done in the exact same way
	* TwinsUK working though, Need to update Turnbaugh
	
* Need to remove MetaHit data and update the respective Figures
	* Done
	
* Need to update manuscript and .R file so that all calculations are made automatically

* Need to add Power calculation to the overall dataset for Mann-Whitney Rank Sum
	* Use R package samplesize
		* n.wilcox.ord(power=, alpha=, t=, p=c(x, y), q=c(x, y))
			* Ralph Scherer (maintainer)
			* Based on Zhao YD, Rahardja D, Qu YongMing. Sample size calculation for the 
			Wilcoxon-Mann-Whitney test adjusting for ties. Statistics in Medicine 2008; 27:462-468
	* Use R custom script to get power
		* Website can be found [here](http://r.789695.n4.nabble.com/How-to-compute-the-power-of-a-wilcoxon-test-td3815616.html)
			* pval <- replicate(1000, wilcox.test(rnorm(366,.0032,.012), rnorm(366,.00042,.0016))$p.value)
			* summary(pval) 
			* sum(pval < .05) (e.g. 979 would mean a power of 97.9%)
	* Use R package pwr
	
* Need to update SequenceProcessing.Rmd since I changed a few things to make data sets more similar
	

