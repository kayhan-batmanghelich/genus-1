#!/bin/bash

results=(XBA XBC XBCCR XCC XC XCCR)
for idx in {0..5}; do
    python plot_cog.py lg_${results[idx]}.pkl ${results[idx]} ${results[idx]}
done
