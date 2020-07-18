# Speckle aberration recovery

Note AW 7/16/2020: by popular demand we are improving our scripts and translating them into python. Stay tuned!


## Introduction

This repository is the open science addition to the research paper "[Extreme ultraviolet microscope characterization using photomask surface roughness](https://www.nature.com/articles/s41598-020-68588-w)" (Gunjala, 2020), where we demonstrate the use of natural roughness to recover the aberration of a imaging system, in that case the SHARP EUV microscope ([sharp.lbl.gov](sharp.lbl.gov))

![alt text](https://raw.githubusercontent.com/gautamgunjala/speckleAberrationRecovery/master/assets/speckle_tf.gif "speckle through focus")


## Content
This repository has four main folders:

* `data`, where data is available, would want to try your own luck
* `matlab`, where the core simulations and experimental results are available
* `python`, where an open-source version of the experimental processing is available.


## Usage

### MATLAB
Simulations of aberration recovery across many randomly-initialized wavefront error functions of a specified rms magnitude can be performed using the script `run.m`. This script was used to generate Figure 6 of our recent paper.

### Python
There is a jupyter notebook.
you can also run a script.

## Examples

## Additional references
The main concept of this paper is based on preliminary research described in the research paper "[Optical transfer function characterization using a weak diffuser](https://doi.org/10.1117/12.2213271)" (Gunjala, 2016) and proof-of-concept experimentation described in the research paper "[Aberration recovery by imaging a weak diffuser](https://doi.org/10.1364/OE.26.021054)" (Gunjala, 2018).


_Gautam Gunjala and Antoine Wojdyla, July 2020_
