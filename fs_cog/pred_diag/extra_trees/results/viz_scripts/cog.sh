#!/bin/bash

results=(XCC XCCR XC XBCCR XBC XBA)
for idx in {0..5}; do
    python plot_cog.py et_${results[idx]}.pkl ${results[idx]} cog_${results[idx]}
done

