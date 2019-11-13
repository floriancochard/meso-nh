# Description
This is a personal install of Meso-NH, a non-hydrostatic model used for research. It is developed by the Laboratoire d'Aérologie (UMR 5560 UPS/CNRS) and by CNRM-GAME (UMR 3589 CNRS/Météo-France). The holders of Meso-NH software are Centre National de la Recherche Scientifique CNRS, Météo-France and Université Paul Sabatier.

When I tried to install Meso-NH, I found the [online tutorial](http://mesonh.aero.obs-mip.fr/mesonh/dir_doc/book3_clean/node3.html) very technical and hard to follow. So I wrote a simple version for Ubuntu distribution available on Medium [at this link](https://medium.com/@floriancochard/install-meso-nh-on-ubuntu-18-04-16f811de5b8). Be aware that cloning this repo is not enough to make it works, you need to follow the tutorial to build the project from source (`make`, `make install`, etc.). The latest version of Meso-NH can be downloaded [at this link](http://mesonh.aero.obs-mip.fr/mesonh54).

# Side note
Don't forget to load a *profile* each time you start a new session (or configure your `.bash_aliases` accordingly). If you haven't created a default profile yet, do it:
```
cd MNH-V5-4-2/src
./configure
. ../conf/profile_mesonh
```
To create a profile with specific settings:
```
export ARCH=LXgfortran        #LXgfortran,LXifort,LXpgi,AIX64,SX8,BGQ (e.g., use Intel "gfortran" compiler on LX=linux Plateform)
export VER_MPI=MPIAUTO     #MPIVIDE,MPIAUTO,MPIICE,MPIINTEL (e.g., use MPI with compiler wrapper 'auto', for computer having this wrapper installed)
export OPTLEVEL=DEBUG         #OPTLEVEL=DEBUG,O2 (e.g., compile in O2, 4 times faster than DEBUG, but less error checks)
./configure`
. ../conf/.profile_mesonh-LXgfortran-R8I4-MNH-V5-4-2-MPIAUTO-DEBUG`
```
