#!/bin/bash

results=(XBA XBCCR XBCOV XBC XBCR XB)
for idx in {0..5}; do
    python plot_brain.py lg_${results[idx]}.pkl ${results[idx]} ${results[idx]}
done

