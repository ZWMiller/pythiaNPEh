#!/bin/csh

source /star/u/huangbc/.cshrc
cd /star/u/zbtang/myTools/root/
source bin/thisroot.csh
cd /star/u/huangbc/data01/simu/pythia8/zbtangBtojpsi
./pmainBJpsiHcorr cards/BJpsi_21.cmnd output/BJpsiHcorr_21.root B
