README file for Matlab source code supporting the paper "Protein Complexes Discovery Based on Protein-Protein Interaction Data via a Regularized Sparse Generative Network Model".

Contents of this archive
------------------------
This archive contains several Matlab scripts used for detecting protein complexes in PPI networks using the 
algorithm RSGNM described in the above paper. 

1) RSGNM.m: Matlab script for the main algorithm of RSGNM described in the Algorithm 1 of the above paper.

2) multi_RSGNM.m: Matlab script that repeats the entire calculation of RSGNM multiple times and choose the result
that gives the lowest value of the obejective function of (8). 

3) demo_RSGNM.m: A simple Matlab script to test RSGNM. When a PPI network is choosed, it can be run in a straightforward 
manner within a Matlab window.


This archive also contains a folder named as "data" which includes the three PPI networks used in the study. The three PPI 
networks are saved with Matlab .mat format. For each network, it is save as a structure that contains information of PPI 
network. Filed of "adjacent matrix" is the adjacent matrix PPI network and filed of "protein_list" is the name list of 
correspding proteins.

Please do not hesitate to contact Dao-Qing Dai at stsddq@mail.sysu.edu.cn (or Xiao-Fei Zhang at zhangfx9@mail2.sysu.edu.cn)
to seek any clarifications regarding any contents or operation of the archive.
