# Quality Threshold Clustering of Molecular Dynamics

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)

Clustering Molecular Dynamics trajectories is a common analysis that allows grouping together similar conformations. Several algorithms have been designed and optimized to perform this routine task and among them, Quality Threshold (QT) stands as a very attractive option. QT guarantees that in retrieved clusters, no pair of frames will have a similarity value greater than a specified threshold and hence a set of strongly correlated frames is obtained for each cluster. For more information about QT, please refer to Heyer et. al. work (Heyer, L. J.; Kruglyak, S.; Yooseph, S. Exploring Expression Data: Identification and Analysis of Coexpressed Genes. 1999, No. 213, 1106–1115. Genome Research).

## Prerequisites

To successfully run QT.py script, you should have previously installed the following dependencies:

* **MDTraj**: Fast pairwise RMSD computation
* **Numpy**: Numerical calculations

## Installing dependencies

* **MDTraj**

You can install MDTraj via pip or conda. The latter is recommended.

``
$ conda install -c conda-forge mdtraj
``

or

``
$ pip install mdtraj
`` 

* **Numpy**

If Numpy is not already present in your system, see the [recommended options](https://scipy.org/install.html) to install it. 


## Getting Started

After installation of dependencies, you can access the help of this script as **python QT.py -h**. QT.py can deal with a variety of topology/trajectory formats and a rich selection syntax is also available. Both features are provided by the MDTraj library, please refer to its [documentation](http://mdtraj.org/1.9.3/) for further details. 

The following line exemplifies available arguments. Note that only the **-traj** is mandatory (if trajectory does not contain topological information, then **-top** argument is also required).

```
python QT.py -traj trajectory.dcd -top topology.psf -first 0 -last 100 -stride 5 -sel backbone -cutoff 4 -minsize 50 -odir ./
```

**-traj**: the path to the trajectory file

**-top**: the path to the topology file

**-first**: first frame to analyze (zero-based counting)

**-last**: last frame to analyze (zero-based counting)

**-stride**: stride (step) of frames to analyze 

**-sel**: selection syntax whose atoms will be superposed before RMSD calculation

**-cutoff**: threshold of similarity (in Angstroms)

**-minsize**: minimum number of frames inside returned clusters

**-odir**: the directory where to save output files

QT.py produces two files: 1) **QT_Clusters.txt** where each line represents a frame and contains the number of the cluster that frame belongs to and 2) **QT_Visualization.log** that is in NMRcluster format and allows users to interface with VMD program.


## VMD interface

If you have VMD program and the [clustering plugin](https://github.com/luisico/clustering) by Luis Gracia, you can
use it for importing the **QT_Visualization.log** file and visualize the results of QT.py (once you have load the trajectory in VMD).

This option is under revision ...

## Authors 
[![Ask Me Anything !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](https://www.linkedin.com/in/roy-gonz%C3%A1lez-alem%C3%A1n-aba3aba3/)

[**Roy González-Alemán**](https://www.linkedin.com/in/roy-gonz%C3%A1lez-alem%C3%A1n-aba3aba3/) 

## License

This project is licensed under the GNU GENERAL PUBLIC LICENSE Version 3

## Citation

[![DOI:10.1021/acs.jcim.9b00558](https://zenodo.org/badge/DOI/10.1021/acs.jcim.9b00558.svg)](https://doi.org/10.1021/acs.jcim.9b00558)

If this work has been useful for your research, please cite it:

González-Alemán, R.; Hernández-Castillo, D.; Caballero, J.; Montero-Cabrera, L. A. Quality Threshold Clustering of Molecular Dynamics: A Word of Caution. J. Chem. Inf. Model. 2019, acs.jcim.9b00558.

