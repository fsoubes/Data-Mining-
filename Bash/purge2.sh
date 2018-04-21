#!/bin/bash

cat cluster.txt | tr -d '(u)'|tr -d \]  | tr -d \[ | tr -d \' > clustered.txt
