#!/usr/bin/perl

while(<>){
    @field=split(/:/,$_);
    system("mv $field[0] tmp");
}

