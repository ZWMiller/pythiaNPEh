#!/bin/csh

source /star/u/zamiller/.cshrc
cd /star/u/zbtang/myTools/root/
source bin/thisroot.csh
cd /star/u/zamiller/simu/NPETemplates
./NPEHDelPhiCorr cards/NpeC_9.cmnd output/NpeCHcorr_9.root C
