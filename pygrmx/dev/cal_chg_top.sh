#!/bin/bash



awk '/atoms/{L=1} L{sum+=$7} /bonds/{L=0} END {printf "total charge : %10.6f\n", sum} ' $1


