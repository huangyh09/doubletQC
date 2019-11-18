#!/bin/sh

dat_dir=../data
bin_dir=../bin

for sample in demux_ctrl demux_stim
do
    python $bin_dir/run_scrublet.py -i $dat_dir/$sample
    Rscript $bin_dir/run_doubletFinder.R $dat_dir/$sample # $dat_dir/$sample/doubletFinder_table.tsv
    
    python $bin_dir/run_combineResults.py \
        --dbFinder $dat_dir/$sample/doubletFinder_table.tsv \
        --scrublet $dat_dir/$sample/scrublet_table.tsv
done