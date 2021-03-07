COMPAS_HPC
------------
Here are some basic instructions for setting up and running COMPAS_HPC. COMPAS is our binary stellar evolution code. COMPAS_HPC is a suite of tools written in python to enable running COMPAS on High Performance Computers (HPC), also known as supercomputers (see supported below). 

Supported HPC Facilities
-------------------------
`OzSTAR <https://supercomputing.swin.edu.au>`_ (Swinburne)

Tsunami (Birmingham)

Requirements
--------------
I will assume here a basic knowledge of linux commands and access to one of the above HPC facilities.

Instructions
---------------
COMPAS is hosted on gitlab here:
http://gitlab.sr.bham.ac.uk/COMPAS/COMPAS
 
1) Sign up for a gitlab account:

You should be able to sign up for a gitlab account; it is best to use a gmail account if you can, other accounts sometimes get denied.
 
2) Add an ssh key (optional)

We need to generate an ssh key on OzSTAR (and possibly also your local machine) to clone the repository. Do this by:
 
`ssh-keygen -t rsa -C "your.email@example.com" -b 4096`

where you replace your email and respond to the prompt. Then, copy the public part of your key using e.g.

`vim ~/.ssh/id_rsa.pub`

and then paste this into a new entry on the gitlab website:

http://gitlab.sr.bham.ac.uk/profile/keys

3) Source COMPAS_SDK

COMPAS relies on a number of packages and libraries. You may want to set up your own environment eventually, but for now on OzSTAR you can source mine:
 
`source /fred/oz101/sstevens/COMPAS_SDK/COMPAS_SDK`

Remember to do this every time you want to run COMPAS.
 
4) Download COMPAS

The next step is to download COMPAS from gitlab. If you made an ssh key in step 2, you can do this by doing:

`git clone git@gitlab.sr.bham.ac.uk:COMPAS/COMPAS.git`

Otherwise, you can do this by using https by doing:

`git clone http://gitlab.sr.bham.ac.uk/COMPAS/COMPAS.git`

and entering your gitlab credentials when requested.
 
5) Set `COMPAS_ROOT_DIR`
 
COMPAS_HPC refers to things using paths relative to the $COMPAS_ROOT_DIR. You can set this by running:
 
`export COMPAS_ROOT_DIR=/path/to/COMPAS`
 
where `/path/to/COMPAS` is the path to the directory you just downloaded COMPAS to.

You need to do this every time you run COMPAS. You may want to add this command to your `~/.bash_profile`.
 
Make sure you source your bash profile by doing:
 
`source ~/.bash_profile`
 
I will now use `$COMPAS_ROOT_DIR` to refer to this directory. This variable is also used by some of the codes.
 
6) Compile COMPAS
 
COMPAS is written mostly in C++ and thus requires compiling. Go to the source directory:

`cd $COMPAS_ROOT_DIR/COMPAS`
 
and compile by doing:
 
make --makefile=Makefile_G2 COMPAS
 
where you use the relevant makefile for the HPC facility you are on. Note that the same makefile is used for OzSTAR as G2. This should compile COMPAS with no errors.
 
7) Test
 
Once finished, test COMPAS by doing e.g.:
 
`./COMPAS --help`
 
8) compasHPC
 
compasHPC is a suite of python tools written to make running COMPAS on High Performance Computers (HPC) very easy. The code is in:
 
`$COMPAS_ROOT_DIR/CompasHPC/`
 
You typically need to edit 2 files:
 
`$COMPAS_ROOT_DIR/CompasHPC/compas_hpc_input.py`
`$COMPAS_ROOT_DIR/CompasHPC/masterFolder/pythonSubmit.py`
 
9) `compas_hpc_input.py`

This contains 'meta settings' such as where you want to send the output of your jobs and how many jobs to run, along with other settings. Make sure to set
 
`cluster = 'ozstar' # or tsunami`
 
On OzSTAR, you need to send your output to fred e.g.
 
`rootOutputDir = '/fred/oz101/sstevens/popsynth/test'`
 
and make sure nBatches is small to start with (say, 3). Because COMPAS is embarassingly parallel, we can gain a factor nBatches speed up by splitting any job in to nBatches pieces, which can be very powerful if nBatches is say 100. It can also cause lots of problems very quickly, so be careful!
 
10) `pythonSubmit.py`

This contains the 'physics settings' that are sent to COMPAS when you run a job. You can choose how many binaries to evolve, what metallicity to use, what mass distribution and so on. There's a lot there, I won't try to explain it all here. Try 

`$COMPAS_ROOT_DIR/COMPAS/COMPAS --help`
 
11) Run `compas_hpc.py`

Finally, run
 
`python $COMPAS_ROOT_DIR/CompasHPC/compas_hpc.py`
 
to run COMPAS. This will automatically submit the job to the relevant job scheduler for your HPC facility. On OzSTAR that is Slurm, on G2 that is PBS and on tsunami that is Condor. There should be some pretty verbose output showing you the commands the code is running.