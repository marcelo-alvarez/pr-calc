#!/bin/bash

scinetgpc_desc="GPC Cluster at Scinet"
scinetgpc_file="scinet-gpc"

scinetbgq_desc="Blue Gene\/Q at Scinet"
scinetbgq_file="scinet-bgq"

cita_desc="CITA Workstations"
cita_file="cita"

sunnyvale_desc="Sunnyvale Cluster at CITA"
sunnyvale_file="sunnyvale"

darwin_desc="Generic Mac OS machine"
darwin_file="darwin"

stampede_desc="Stampede Cluster at TACC"
stampede_file="tacc-stampede"

for f in scinetgpc scinetbgq cita sunnyvale darwin stampede
do

    desc_var=`echo $f\_desc`
    file_var=`echo $f\_file`

    desc=${!desc_var}
    file=${!file_var}

    fout=Settings.$file

    cat machines/directions machines/$file > $fout

    sed -i "s/DESC_REPLACE/$desc/g" $fout
    sed -i "s/FILE_REPLACE/$fout/g" $fout

done


