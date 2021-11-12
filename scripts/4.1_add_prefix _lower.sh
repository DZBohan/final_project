#!/bin/sh
for files in $(ls *.count)
    do mv $files "lowerlobe-"$files
done
