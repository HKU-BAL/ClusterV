import pysam
import argparse
import statistics
import sys
import os


def main(args):
    alm_list=[]
    with pysam.AlignmentFile(args.bam, "rb" ) as samfile, pysam.AlignmentFile(args.out_fn+'.tagged.bam','wb',template=samfile) as fo:
        for read in samfile.fetch():
            read.set_tag("X?",args.tag,"Z")
            alm_list.append(read)
        for alm in alm_list:
            fo.write(alm)
    os.system(f'samtools index {args.out_fn}.tagged.bam')


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', metavar='bam_fn', type=str,)
    parser.add_argument('tag', metavar='tag', type=str,)
    parser.add_argument('--out_fn', type=str)

    args = parser.parse_args()
    main(args)
