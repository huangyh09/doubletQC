doubletQC
=========

Benchmarking of doublet detection methods in single-cell RNA-seq data


Method list
-----------
* `Scrublet <https://github.com/AllonKleinLab/scrublet>`_

  **Model parameter:**
  
  * r: default 2. The ratio of simulated cells
  * K: sqrt(n_input_cell) / 2. The number of neighbours in KNN
  
  **Threshold parameter:**
  
  * doublet score: threshold to set. Can be optimised by Bayesian Gaussian 
    mixture model, added in the wrapper function here.

* `DoubletFinder <https://github.com/chris-mcginnis-ucsf/DoubletFinder>`_

  **Model parameter:**
  
  * pN: default 0.25. The ratio of simulated cells r = pN / (1 - pN). This 
    parameter pN has been shown resistant in the paper
  * pK: propotion of neighbours in KNN. Can be optimised with build-in function
  
  **Threshold parameter:**
  
  * nExp: number of expected to *heterotypic* doublets for threshold
  * pANN: fraction of simulated doublet neighbours, threshold with nExp

* `DoubletDecon <https://github.com/EDePasquale/DoubletDecon>`_


Wrapper functions
-----------------
See wrapper functions in the `bin folder
<https://github.com/huangyh09/doubletQC/tree/master/bin>`_


Example
-------
* data from Demuxlet paper (`Kang et al, 2018, Nature Biotech
  <https://www.nature.com/articles/nbt.4042>`_)
* bash script: `run_demuxlet_data.sh 
  <https://github.com/huangyh09/doubletQC/blob/master/examples/run_demuxlet_data.sh>`_
* jupyter notebook for analysis: `demuxlet_dataset.ipynb
  <https://github.com/huangyh09/doubletQC/blob/master/examples/demuxlet_dataset.ipynb>`_