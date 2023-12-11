import sys
import os
import shlex
from subprocess import PIPE
from argparse import ArgumentParser
from shared.utils import subprocess_popen, vcf_candidates_from, _run_command
from shared.interval_tree import bed_tree_from, is_region_in
import subprocess
from cv.filter_large_indel import filter_read


# filtering original bam to ensure read cover both end of target region
# run clair ensemble to get initial variants list
def initial_run(args):
    _bam = args.bam_fn
    _bed = args.bed_fn
    _ref = args.ref_fn
    _sample_id = args.sample_id
    _out_dir = args.out_dir
    _cn_threads = args.clair_ensemble_threads

    cmd = 'mkdir -p %s' % (_out_dir)
    _run_command(cmd)

    cmd = 'mkdir -p %s/%s.v' % (_out_dir, _sample_id)
    _run_command(cmd)
    tree = bed_tree_from(bed_file_path=_bed)
    is_tree_empty = len(tree.keys()) == 0
    # read contig name
    _contig = [k for k, v in tree.items()][0]
    _l, _r = [(list(v)[0][0], list(v)[0][1]) for k, v in tree.items()][0]
    _end_dis, _end_len = 400, 200
    _l_range = _l + _end_dis, _l + (_end_dis + _end_len)
    _r_range = _r - (_end_dis + _end_len), _r - _end_dis
    _contig, _l_range, _r_range

    # get read cover bed region
    cmd = "samtools view %s -o %s/%s_a.bam -b %s:%s-%s; samtools index %s/%s_a.bam; samtools view %s/%s_a.bam -o %s/%s_o.bam -b %s:%s-%s; samtools index %s/%s_o.bam" % \
        (_bam, _out_dir, _sample_id, _contig, _l_range[0], _l_range[1], _out_dir, _sample_id, _out_dir, _sample_id, _out_dir, _sample_id, _contig, _r_range[0], _r_range[1], _out_dir, _sample_id)
    _run_command(cmd)

    # rm tmp file
    cmd = "rm %s/%s_a.bam" % (_out_dir, _sample_id)
    _run_command(cmd, False)
    cmd = "rm %s/%s_a.bam.bai" % (_out_dir, _sample_id)
    _run_command(cmd, False)

    indel_l = args.indel_l
    if int(indel_l) > 0:
        input_f = "%s/%s_o.bam" % (_out_dir, _sample_id)
        out_f = "%s/%s_f.bam" % (_out_dir, _sample_id)
        filter_read(input_f,  out_f, indel_l)


        cmd = "rm %s/%s_o.bam" % (_out_dir, _sample_id)
        _run_command(cmd, False)
        cmd = "rm %s/%s_o.bam.bai" % (_out_dir, _sample_id)
        _run_command(cmd, False)
    else:
        cmd = "mv %s/%s_o.bam %s/%s_f.bam" % (_out_dir, _sample_id, _out_dir, _sample_id)
        _run_command(cmd, False)
        cmd = "mv %s/%s_o.bam.bai %s/%s_f.bam.bai" % (_out_dir, _sample_id, _out_dir, _sample_id)
        _run_command(cmd, False)


    print('running clair ensenmble')
    _py_s_d = os.path.dirname(os.path.abspath(__file__))
    run_clair_path = "%s/run_Clair_ensemble_cv.sh" % (_py_s_d)

    cmd = "bash %s %s/%s_f.bam %s %s %s %s/%s.v %s > %s/%s.v/run.log" % \
    (run_clair_path, _out_dir, _sample_id, _sample_id, _ref, _bed, _out_dir, _sample_id, _cn_threads, _out_dir, _sample_id)
    _run_command(cmd)


def main():
    parser = ArgumentParser(description="Initialize ClusterV setting")

    parser.add_argument('--out_dir', type=str, default="output_dir",
                    help="Output folder, required")

    parser.add_argument('--bam_fn', type=str, default="input.bam",
                        help="input bam file, required")

    parser.add_argument('--bed_fn', type=str, default="input.bed",
                        help="input target regions bed file, required")

    parser.add_argument('--ref_fn', type=str, default="ref.bed",
                        help="input reference fasta, required")

    parser.add_argument('--sample_id', type=str, default="hiv_sample_1",
                        help="sample_id, optional")

    parser.add_argument('--indel_l', type=int, default=50,
                        help="filtering read with indel length > indel_l [50], set [0] to disable filtering, optional")

    parser.add_argument('--clair_ensemble_threads', type=int, default=16,
                        help="Clair-Ensemble threads, we recommend using 16, [16] optional")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    initial_run(args)

if __name__ == "__main__":
    main()
