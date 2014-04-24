 #!/bin/bash
headd=$1
NjobS=$2
NjobE=$3
NjobCur=$NjobS

while [ $NjobCur -le $NjobE ]
do 

    sed -e "s/REPL/$NjobCur/"  < ${headd}_template.qsub > ${headd}_$NjobCur.qsub
    qsub ${headd}_$NjobCur.qsub
    NjobCur=$[$NjobCur+1]

done 
