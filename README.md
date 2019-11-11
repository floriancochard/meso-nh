# Quick start
Cloning the `Meso-NH` repo is not enough to make it works, you need to install dependencies as well as OpenMPI. I've written a short tutorial on Medium to help you install Meso-NH on Ubuntu:
https://medium.com/@floriancochard/install-meso-nh-on-ubuntu-18-04-16f811de5b8

Once you've installed it, don't forget to load a `profile` each time you start a new session (or configure your `.bash_aliases` accordingly). If you haven't already created a profile, do so:
```
cd MNH-V5-4-2/src
./configure
. ../conf/profile_mesonh
```
If you want to set specific compiler, MPI or optlevel:
```
export ARCH=LXgfortran        #LXgfortran,LXifort,LXpgi,AIX64,SX8,BGQ (e.g., use Intel "gfortran" compiler on LX=linux Plateform)
export VER_MPI=MPIAUTO     #MPIVIDE,MPIAUTO,MPIICE,MPIINTEL (e.g., use MPI with compiler wrapper 'auto', for computer having this wrapper installed)
export OPTLEVEL=DEBUG         #OPTLEVEL=DEBUG,O2 (e.g., compile in O2, 4 times faster than DEBUG, but less error checks)
./configure`
```
This will create a new profile in `/conf/`, load it with:
`.profile_mesonh-LXgfortran-R8I4-MNH-V5-4-2-MPIAUTO-DEBUG`
