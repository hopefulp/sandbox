#!/bin/bash

for x in `ls *.log`; do echo $x; grep "F=" $x | tail -5; done

