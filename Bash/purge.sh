#!/bin/bash

cat cluster.txt | tr -d '(u)'|tr \] \\n | tr -d \' | sed -r 's/[/[]+/cluster:  /'  > clusterclean.txt
