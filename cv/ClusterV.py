import sys
import os
import shlex
from subprocess import PIPE
from argparse import ArgumentParser
from shared.utils import subprocess_popen, vcf_candidates_from, _run_command
from shared.interval_tree import bed_tree_from, is_region_in
import subprocess
from cv.run_initial_call import initial_run
from cv.run_ClusterV import CV_run
from cv.get_consensus import run_get_consensus

def main():
    parser = ArgumentParser(description="run ClusterV pipeline")

    # general arguments 
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

    parser.add_argument('--threads', type=int, default=48,
                        help="running threads, we recommend using 48 or above, [48] optional")

    parser.add_argument('--clair_ensemble_threads', type=int, default=16,
                        help="Clair-Ensemble threads, we recommend using 16, [16] optional")

    parser.add_argument('--subtype_parallel', type=int, default=3,
                        help="[EXPERIMENTAL] number of sutypes parallel run Clair, [3] optional")

    # initial filtering options
    parser.add_argument('--indel_l', type=int, default=50,
                        help="filtering read with indel length > indel_l [50], set [0] to disable filtering, optional")

    # clusterV options
    parser.add_argument('--top_k', type=int, default=25,
                        help="top k sutypes to output, optional")

    parser.add_argument('--n_min_supports', type=int, default=50,
                        help="minimum read support for creating a subtype, optional")

    parser.add_argument('--n_max_candidates', type=int, default=15,
                        help="[EXPERIMENTAL] number of selected candidates for clustering, optional") 

    parser.add_argument('--min_af', type=float, default=0.05,
                        help="[EXPERIMENTAL] minimum AF when cluastering, optional")

    parser.add_argument('--n_max_coverage', type=int, default=10000,
                        help="[EXPERIMENTAL] max read for clustering, optional")

    # generate consensus and get HIVDB report options
    parser.add_argument('--hivdb_url', type=str, default="",
                        help="hivdb url defalut query from internet, for localize the HIVDB, please check https://github.com/hivdb/sierra, and update this setting accordingly, e.g. \
                        by using --hivdb_url http://localhost:8111/sierra/rest/graphql")

    parser.add_argument('--number_of_read_for_consense', type=int, default=1000,
                        help="[EXPERIMENTAL] number of original read for generating consense")

    parser.add_argument('--flye_genome_size', type=str, default="5k",
                        help="[EXPERIMEANTAL], flye --genome-size for generating consensus, we recommand using 5k for HIV genome")
    
    parser.add_argument('--flye_genome_size_olp', type=str, default="1000",
                        help="[EXPERIMEANTAL], flye -m for generating consensus, we recommand using 1000 for HIV genome")

    parser.add_argument('--flye_nano_type', type=str, default="nano-hq",
            help="[EXPERIMEANTAL], flye option for different ont type, default --nano-hq, check https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md")

    args = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    if args.subtype_parallel < 1:
        args.subtype_parallel = 1
    if args.threads < 8:
        args.threads = 8

    _out_dir = args.out_dir
    _sample_id = args.sample_id
    _threads = args.threads
    _subtype_parallel = args.subtype_parallel
    args.clair_ensemble_threads = int(args.threads / args.subtype_parallel)

    # run initial run
    print("-----------------------------\n[ ** STEP 1 ** ] begin initial run\n")
    initial_run(args)

    # run clusterV
    args.vcf_fn = "%s/%s.v/out.vcf" % (_out_dir, _sample_id)
    args.bam_fn = "%s/%s_f.bam" % (_out_dir, _sample_id)
    args.out_dir = "%s/clustering" % (_out_dir)
    print("-----------------------------")
    print("[ ** STEP 2 ** ] begin clusterV run, with [%s threads, parallel %s subtypes with %s clair-ensemble threads]\n" % (args.threads, args.subtype_parallel, args.clair_ensemble_threads))
    CV_run(args)

    # get consensus and get HIVDB report
    print("-----------------------------\n[ ** STEP 3 ** ]get consensus and HIVDB report\n")
    args.out_dir = "%s/consensus" % (_out_dir)
    args.tar_tsv = "%s/clustering/all_clusters_info.tsv" % (_out_dir)
    run_get_consensus(args)


if __name__ == "__main__":
    main()
