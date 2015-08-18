#!/bin/csh

source /star/u/zamiller/.cshrc
cd /star/u/zbtang/myTools/root/
source bin/thisroot.csh
cd /star/u/zamiller/simu/NPETemplates
./NPEHDelPhiCorr cards/NpeC_28.cmnd output/NpeCHcorr_28.root C
