#!/bin/bash
#  Glenn Gaffield - April 2020
#
#  Fix sequence coordinates in a .tblout file

awk '{if (/^#/) {print} else {
        match($1,"([^/]+)/([0-9]+)-([0-9]+)",accno);          # parse target accession number
        start= accno[2]; end= accno[3]; accno_len= RLENGTH;   # accno[1]/accno[2]-accno[3]
        match($0,"( +"$8")( +"$9")( +\\"$10")",pos);          # save "from", "to", and "strand" fields
        if (start<end) {from=start+$8-1; to=start+$9-1}       # target has forward notation
        else {from=start-$9+1; to=start-$8+1};              # target has reverse notation
        new_from = sprintf("%"length(pos[1])"d",from);        # make new "seq from" field
        new_to = sprintf("%"length(pos[1])"d",to);            # make new "seq to" field
        sub("+","\\+",pos[0]);                                # convert seq_from+seq_to+strand into regexp
        sub(pos[0],new_from""new_to""pos[3],$0);              # string replace seq_from+seq_to+strand with new
        new_accno=sprintf("%-"accno_len"s",accno[1]);         # make new target accession without coords
        sub(accno[0],new_accno,$0);                           # string replace target accno with new
        print;
    }
}' $1
