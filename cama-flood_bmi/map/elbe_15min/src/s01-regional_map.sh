#!/bin/sh

# Please edit "region_info.txt" to make a regional map

mkdir ../hires

./cut_domain
./combine_hires
./generate_inpmat-1deg
./set_map


