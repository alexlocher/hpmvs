# Hierarchical Progressive Multiview Stereo

**Author:** [Alex Locher](http://www.vision.ee.ethz.ch/~alocher)

Progressive multiview stero is an implementation of an algorihm taking a sparse 
3D model (result of structure from motion software) and outputs a dense set of 
surface patches. Due to its progressiveness, the output gets more and more accurate 
the longer it runs.

#### Multiview Dense Stereo
<img src="http://www.vision.ee.ethz.ch/~alocher/pdf/cvpr16/mvs_problem.jpg" 
alt="Tsukuba Dataset" width="80%" border="10" />

#### Progressive Output
<img src="http://www.vision.ee.ethz.ch/~alocher/pdf/cvpr16/progressive_mvs.png" 
alt="Tsukuba Dataset" width="80%" border="10" />


## Related Publication
[1] Alex Locher, Michal Perdoch and Luc Van Gool. Progressive prioritized multi-view stereo. *CVPR 2016*. 

[2] Y. Furukawa and J. Ponce, “Accurate, dense, and robust multiview stereopsis,” IEEE Trans. Pattern Anal. Mach. Intell., vol. 32, no. 8, pp. 1362–1376, 2010.


# Licence
HPMVS is realeased under the [GPLv3 Licence](https://www.gnu.org/licenses/gpl-3.0.txt). 


If you use the algortihm in your academic work, please cite:

    @article{locher16,
      title={Progressive Prioritized Multi-view Stereo},
      author={Locher, Alex and Perdoch, Michal and Van Gool, Luc},
      journal={CVPR},
      year={2016}
     }


# Build and Install
We use the cmake framework for compilation and installation. The software is currently
only tested in a Linux environment. 

```
git clone https://github.com/alexlocher/hpmvs
cd hpmvs
mkdir build && cd build
cmake ..
make -j4
```

# Usage
The binary *hpmvs* takes an [nvm-file](http://ccwu.me/vsfm/doc.html#nvm) as input and you have to specify an output directory:

```
./hpmvs --nvm=<nvm-file> --outdir=/tmp
```

