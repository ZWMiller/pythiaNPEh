#!/bin/csh

source /star/u/huangbc/.cshrc
cd /star/u/zbtang/myTools/root/
source bin/thisroot.csh
cd /star/u/huangbc/data01/simu/pythia8/zbtangBtojpsi
./pmainBJpsiHcorr cards/BJpsi_37.cmnd output/BJpsiHcorr_37.root B
