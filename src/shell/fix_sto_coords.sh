#!/bin/bash
#  Glenn Gaffield - April 2020
#
#  Fix sequence coordinates in a .sto file

awk '{re="([^ \t,;=]+)/([0-9]+)-([0-9]+)/([0-9]+)-([0-9]+)([= ]+)"; 
    while (match($0,re,accno)>0) {
      if ($0~"^#=GF DUPL") {RLENGTH=""};
      start=accno[2]; end=accno[3];
      from=accno[4]; to=accno[5];
      if (start<end) {new_start=start+from-1; new_end=start+to-1}       # accno has forward notation
        else {new_start=start-from+1; new_end=start-to+1};              # accno has reverse notation
      new_accno=sprintf("%-"RLENGTH"s",accno[1]"/"new_start"-"new_end""accno[6]);
      sub("\\|","\\|",accno[0]);                                        # escape the "|" character
      sub(accno[0],new_accno,$0);
    };
    print; }' $1
