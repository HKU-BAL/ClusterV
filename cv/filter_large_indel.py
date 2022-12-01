import argparse
import sys
import os
import pysam
import re
def filter_read(input_f, out_f, indel_l = 50):
    print("input: %s\noutput: %s\nfiltering indel length: %s" % (input_f, out_f, indel_l))

    with pysam.AlignmentFile(input_f, "rb" ) as samfile, pysam.AlignmentFile(out_f,'wb',template=samfile) as fo:
        aln_list = []
        for read in samfile:
            #print(read.cigarstring)
            tar_ciagr = read.cigarstring
            #for num1, i_or_d, num2, m in re.findall('(\d+)([ID])(\d+)?([A-Za-z])?', tar_ciagr):
            is_has_large_indel = False
            for num1, i_or_d in re.findall('(\d+)([ID])', tar_ciagr):
                if int(num1) > indel_l:
                    is_has_large_indel = True
                    #print(num1, i_or_d)
            if is_has_large_indel:
                #print(tar_ciagr)
                continue
            aln_list.append(read)
        for aln in aln_list:
            fo.write(aln)
    os.system(f'samtools index {out_f}')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_f', type=str, default='bam_fn')
    parser.add_argument('--out_f', type=str, default='out_f')
    parser.add_argument('--indel_l', type=int, default=50)
    args = parser.parse_args()
    input_f = args.input_f
    out_f = args.out_f
    indel_l = args.indel_l
    filter_read(input_f, out_f, indel_l)

if __name__=='__main__':
    main()
