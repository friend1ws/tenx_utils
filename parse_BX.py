#! /usr/bin/env python

import sys, os, pysam

def parse_BX_tag(input_bam, output_file_handle, chr, start, end, additional_str):

    bamfile = pysam.Samfile(input_bam, 'rb')
    hout = output_file_handle
 
    # maybe add the regional extraction of bam files
    for read in bamfile.fetch(chr, int(start) - 1, int(end) - 1):

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip if either of the read pair is unmapped
        if flags[2] == "1" or flags[3] == "1": continue

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        BX_current = None
        for item in read.tags:
            if item[0] == "BX": BX_current = item[1]
                    
        if BX_current is None: continue

        # get the alignment basic information
        chr_current = bamfile.getrname(read.tid)
        pos_current = str(read.pos + 1)
        dir_current = ("-" if flags[4] == "1" else "+")
        read_name = read.qname + ("/1" if flags[6] == "1" else "/2")
        mapq_current = str(read.mapq)

        print >> hout, read_name + '\t' + chr_current + '\t' + pos_current + '\t' + dir_current + '\t' + mapq_current + '\t' + BX_current + '\t' + additional_str



input_sv_list = sys.argv[1]
input_bam = sys.argv[2]
output_file = sys.argv[3]

hout = open(output_file, 'w')
with open(input_sv_list, 'r') as hin:

    header_line = "#"
    while header_line.startswith("#"):
        header_line = hin.readline().rstrip('\n')

    print >> hout, header_line
    header = header_line.split('\t')

    header2ind = {}
    for i in range(0, len(header)):
        header2ind[header[i]] = i

    for line in hin:
        F = line.rstrip('\n').split('\t')
        if float(F[header2ind["Minus_Log_Fisher_P_value"]]) < 2: continue

        sv_key = F[0] + ':' + F[2] + F[1] + '-' + F[3] + ':' + F[5] + F[4] + ',' + F[6]

        if F[2] == "+":
            parse_BX_tag(input_bam, hout, F[0], int(F[1]) - 10000, int(F[1]), sv_key + '\t' + '1')
        else:
            parse_BX_tag(input_bam, hout, F[0], int(F[1]), int(F[1]) + 10000, sv_key + '\t' + '1')

        if F[5] == "+":
            parse_BX_tag(input_bam, hout, F[3], int(F[4]) - 10000, int(F[4]), sv_key + '\t' + '2')
        else:
            parse_BX_tag(input_bam, hout, F[3], int(F[4]), int(F[4]) + 10000, sv_key + '\t' + '2')





