#!/bin/csh

source /star/u/zamiller/.cshrc
cd /star/u/zbtang/myTools/root/
source bin/thisroot.csh
cd /star/u/zamiller/simu/NPETemplates
./NPEHDelPhiCorr cards/NpeB_25.cmnd output/NpeBHcorr_25.root B
