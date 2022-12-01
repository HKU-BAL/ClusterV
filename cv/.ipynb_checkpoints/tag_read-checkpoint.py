import pysam
import argparse
import statistics
import sys
import os
from cv.get_list_from_file import get_dict_from_csv


def run_tag_read(_read_tag_csv, _n_cluster, _bam, _out_fn, _split_bam):
    read_group_dict=get_dict_from_csv(_read_tag_csv)
    alm_list=[]
    if _split_bam:
        #  Xtag.csv
        #n_cluster=int(args.read_tag_csv[args.read_tag_csv.find('tag.csv')-1])
        n_cluster=_n_cluster
        split_alms_dict={ str(i+1): [] for i in range(n_cluster)}

    with pysam.AlignmentFile(_bam, "rb" ) as samfile, pysam.AlignmentFile(_out_fn+'.tagged.bam','wb',template=samfile) as fo:
        for read in samfile.fetch():
            if read.qname not in read_group_dict:
                continue
            read.set_tag("X?",read_group_dict[read.qname],"Z")
            alm_list.append(read)
            if _split_bam:
                split_alms_dict[read_group_dict[read.qname]].append(read)

        for alm in alm_list:
            fo.write(alm)

        if _split_bam:
            for k,v in split_alms_dict.items():
                with pysam.AlignmentFile(f'{_out_fn}.tagged.tag{k}.bam','wb',template=samfile) as split_bam_out:
                    for alm in v:
                        split_bam_out.write(alm)
            for i in range(n_cluster):
                os.system(f'samtools index {_out_fn}.tagged.tag{i+1}.bam')
    os.system(f'samtools index {_out_fn}.tagged.bam')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', type=str, default='bam_fn')
    parser.add_argument('--read_tag_csv', type=str, default='read_tag_csv_fn')
    parser.add_argument('--out_fn', type=str, default='out_fn')
    parser.add_argument('--n_cluster', type=int, default=2)
    parser.add_argument('--split_bam', default=False, action='store_true')
    parser.set_defaults(split_bam=False)

    args = parser.parse_args()
    
    _read_tag_csv = args.read_tag_csv
    _n_cluster = args.n_cluster
    _bam = args.bam
    _out_fn = args.out_fn
    _is_split_bam = True if args.split_bam else False
  
    run_tag_read(_read_tag_csv, _n_cluster, _bam, _out_fn, _is_split_bam)

    
if __name__=='__main__':
    main()
