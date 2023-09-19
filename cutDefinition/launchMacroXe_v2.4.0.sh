#!/bin/bash

while read line; do
  if [[ $line != *"#"* ]] && [[ ! (-z "$line") ]]; then
    echo $line
    ./optimizeCuts_Xe_v2.4.0 $line
  fi
done <runListXe.txt
