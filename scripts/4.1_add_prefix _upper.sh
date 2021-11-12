#!/bin/sh
for files in $(ls *.count)
    do mv $files "upperlobe-"$files
done
