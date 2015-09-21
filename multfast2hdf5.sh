#!/bin/bash

LABEL1='0'
LABEL2='0'
LABEL3='0'
LABEL4='0'

OPTIONS_FILE=''

OUTPUT_FILE=''

function show_help {
echo "Convert multiple .fast files into the HDF5 format."
echo "Options:"
echo "-h: Show help"
echo "-l: Labels for the channels on the used inputs. Usage: l1,l2,l3,l4"
echo "-o: Specify output file"
echo "-p: Specify further options"
}

while getopts ":hl:o:p:" opt; do
    case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    l)  arrIN=(${OPTARG//,/ })
        LABEL1=${arrIN[0]}
        LABEL2=${arrIN[1]}
        LABEL3=${arrIN[2]}
        LABEL4=${arrIN[3]}
        ;;
    o)  OUTPUT_FILE=$OPTARG
        ;;
    p)  OPTIONS_FILE=$OPTARG
        ;;
    :)  echo "Option -$OPTARG requires an argument." >&2
        exit 1
    esac
done

if [ "$OUTPUT_FILE" = "" ]; then
echo "Output file needs to be specified by -o output_file.h5"
exit 0
fi

shift $((OPTIND-1))

for f in $@
do
  echo "File $f"
  echo "Label1: $LABEL1 LABEL2: $LABEL2"
  fast2hdf5 $f $OUTPUT_FILE $LABEL1 $LABEL2
  echo "Conversion of $f completed"
done

