# README

for the model associated with: 

Sanda P, Hlinka J, van den Berg M, Bazhenov M, Keliris GA, Krishnan G: Cholinergic modulation supports dynamic switching of resting state networks through selective DMN suppression (2024). *PLoS Comp 
Bio*, accepted. [doi:10.1371/journal.pcbi.1012099](https://doi.org/10.1371/journal.pcbi.1012099)

This model is written in C++. Tested on Debian 11 (GCC 10) & CentOs 8 (GCC 8), but should work on widely on other platforms.
By default set to use 40 threads (change `NUM_THREADS` in `Kncl.cpp` accordingly to your hardware).

Makefile provides instructions for building and running of resting baseline and ACh-modulated activity (see Fig 4, row B of the original paper).

Long range connectivity - DMN of rat manually extracted from NeuroVIISAS atlas of rodent brain (Schmitt and Eipert 2012).

There are more outputs, but the mainly useful ones are in the following files:
- `spikes` - individual spikes
- `time_cx_new` - voltage of pyramidal neurons
- `time_cx_syn` - synaptic input for all cells
- `time_K` - ionic levels (K, Na, Cl, Ca)
