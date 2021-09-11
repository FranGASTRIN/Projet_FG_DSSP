# Projet_FG_DSSP
Script for secondary structure assignment using hydrogen bonds.

Folder downloadable by entering the following command in your bash terminal:

`git clone https://github.com/FranGASTRIN/Projet_FG_DSSP.git`
## Requirements
### For the user who does not have Miniconda (advisable)
Miniconda can be downloaded by going to the following address: https://docs.conda.io/en/latest/miniconda.html

Then enter the following command in your shell, if you are on Linux: `bash Miniconda3-latest-Linux-x86_64.sh`

And if you are on MacOSX :`bash Miniconda3-latest-MacOSX-x86_64.sh`
### Libraries
For the proper functioning of the program, it is recommended to have the following module installed:
- scipy

If Anaconda or Miniconda is installed on your computer, you can directly create the appropriate environment located in the yaml file of this Github repository by entering the following command:
`conda env create -f projectFG.yml`

You will be able to activate the created environment by entering next:
`conda activate FG_dssp`

If you do not have a conda environment, it is still possible to install this module by entering the command:
`python -m pip install --user scipy`
### Processed PDB files
This program takes input of the PDB files in which the hydrogen atoms have been added. If you have basic PDB files, you can add the hydrogens by going to the [Molprobity](http://molprobity.biochem.duke.edu/index.php) website and following the instructions.

If you have a conda environment, it is possible to install the "reduce" command that will allow you to add the hydrogens using
`conda install -c ostrokach-forge reduce`. The command is used as in the following example:

`reduce -BUILD -i 5jjt.pdb > 5jjtFH.pdb`
## Usage
The FG_dssp program takes 2 arguments as input. The first argument is the name of the PDB file containing the hydrogen atoms. The second is the letter, in capital letters, of the chain concerned. Basically, the program launches as the following way :

`python3 FG_dssp.py 5jjtFH.pdb B`

If your files are in different directories, as is the case in this Github repository, it is still possible to enter the command this way:

`python3 ./Projet_FG_DSSP/src/FG_dssp.py ./Projet_FG_DSSP/data/5jjtFH.pdb B`

The results of the program are displayed by default in the terminal. If you want to keep the results obtained in a file, you can use ">" to transfer them to a new file. 

Example: `python3 FG_dssp.py 5jjtFH.pdb B > 5jjt_dssp.txt`
