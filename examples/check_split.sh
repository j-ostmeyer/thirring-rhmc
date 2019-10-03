#!/bin/bash

for filename in fort.200 fort.100 fort.11 control output
do
    diff <(cat $filename.*) $filename || echo "Problem with the pieces of $filename..."
done
