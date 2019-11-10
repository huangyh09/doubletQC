#!/bin/sh

dat_dir=../data
bin_dir=../bin

for sample in demux_ctrl demux_stim
do
    python $bin_dir/run_scrublet.py -i $dat_dir/$sample
    # Rscript $bin/run_doubletFinder.R $dat_dir/$sample
done