
#! /bin/bash
for i in {0..49}
do
cp bingchuExamples/ccbar_$i.cmnd ./NpeC_$i.cmnd
echo "NpeC_$i.cmnd moved"
done
