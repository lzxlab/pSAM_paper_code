# pSAM is a Python package based on deep learning model to perform precise prediction of protein subcellulcar localization
# Overview
Many shuttling proteins transport between different subcellular locations in response to external signals and regulate a broad spectrum of biological processes . This aberrant localization of proteins is often driven by binding to specific molecules that recognize distinct targeting signals. The dynamic nucleocytoplasmic transport process is regulated by the nuclear localization signal (NLS) and nuclear export signal (NES), while mitochondrial and extracellular localization are determined by pre-sequences such as mitochondrial targeting sequence (MT) and signal peptide (SP). In recent years, an increasing number of proteins initially found in one subcellular component have also been observed to localize in other components and perform completely different functions, especially in tumor and disease states. These discoveries further underscore the importance of further investigations on the regulation of protein subcellular localization. 
# Install and use
pSAM could be installed from GitHub. [conda](https://anaconda.org/anaconda/conda) is required to easily install the package. A webserver version of this model could be accessed from http://inuloc.omicsbio.info/.
```
git clone https://github.com/zhengyq1/pSAMmodel
cd pSAMmodel
conda create -n dl37
conda env create -f environment.yaml
```
Once the pakackage was downloded, the users should modify and write the `path to pSAM` to `predict.py` script.
# System requirements
**Hardware requirements:** `pSAM` package requires a standard computer with enough RAM and with/without GPU.<br>
**Software requirements:** This package is supported for macOS and Linux. The package has been tested on the following systems: CentOS (el8) and macOS (10.14.1).<br>
**Other requirements:** [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit) is also needed if you want to run the model based on NVIDIA GPU.
# Python Dependencies
```
astunparse=1.6.3
attention=4.0
efficientnet=0.0.4
gast=0.3.3
keras-self-attention=0.49.0
numpy=1.18.5
scipy=1.4.1
tensorboard=2.4.1
tensorboard-plugin-wit=1.8.0
tensorflow=2.3.0
tensorflow-estimator=2.3.0
absl-py=0.12.0
aiohttp=3.7.4
astor=0.8.1
async-timeout=3.0.1
blinker=1.4
brotlipy=0.7.0
certifi=2020.12.5
cffi=1.14.5
chardet=3.0.4
coverage=5.5
cryptography=3.4.6
```
# How to run
The pSAM is a easy-to-use command-line package.The model could be run by the following command:
```
cd path_to_pSAM
conda activate dl37
python predict.py input_fasta_file output_path
```
Two parameters are needed: `input_fasta_file` is a common fasta file with protein name started with '>' and protein sequence; `output_path` is the path to work and write and results. The `inputFile.fasta` in the fold is an example.

# Example of input file and output result
The `example_output` fold containing the examples of input file and output results: the `inputFile.fasta` is an example of fasta file pf input; the `merged.data.txt` is the merged output results; the `runInfo.txt` file records the running information.
# Example running code
```
cd path_to_pSAM
conda activate dl37
python predict.py example_output/inputFile.fasta example_output/
```
