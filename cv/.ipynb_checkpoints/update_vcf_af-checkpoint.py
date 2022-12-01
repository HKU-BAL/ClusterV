import pysam
import argparse
import shlex
import statistics
import sys
import os
from shared.interval_tree import bed_tree_from, is_region_in
from shared.utils import subprocess_popen, vcf_candidates_from
from subprocess import PIPE
import subprocess
from collections import Counter, defaultdict

# read vcf
# return a list of [chr, pos, ref_base, alt_base, qual, info_id, info]
def read_vcf(vcf_fn, is_tree_empty=True, tree=None, is_snp_only=False):
    vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))
    vcf_list = []
    v_cnt = {'snp': 0, 'indel': 0}
    for row in vcf_fp.stdout:
        columns = row.strip().split()
        if columns[0][0] == "#":
            continue
        # position in vcf is 1-based
        ctg_name, position = columns[0], columns[1]
        if not (is_tree_empty or is_region_in(tree, ctg_name, int(position))):
            continue
        i = columns[:]
        # chr, pos, ref_base, alt_base, qual, info, info, af
#         tar_info = [i[0], int(i[1]), i[3], i[4], i[5], i[-2], i[-1], float(i[-1].split(':')[-1])]
        # chr, pos, ref_base, alt_base, qual, af
        tar_info = [i[0], int(i[1]), i[3], i[4], i[5], float(i[-1].split(':')[-1])]
        if len(i[3]) == 1 and all([len(_j) == 1 for _j in i[4].split(',')]) == 1:
            v_cnt['snp'] += 1
        else:
            v_cnt['indel'] += 1
        if is_snp_only:
            if len(i[3]) != 1 or len(i[4]) != 1:
                continue
        vcf_list.append(tar_info)
    # sort by chromosome and position
    vcf_list.sort(key=lambda x: (x[0], x[1]))
    return vcf_list, [v_cnt['snp'], v_cnt['indel']]

# take unpaded alignment https://samtools.github.io/hts-specs/SAMv1.pdf
def parse_pileup_out(pileup_bases):
    base_idx = 0
    base_list = []
    _cov = 0
    while base_idx < len(pileup_bases):
        base = pileup_bases[base_idx]
        if base in "ACGTNacgtn#*":
            base_list.append(base)
            _cov += 1
        elif base == '+' or base == '-':
            base_idx += 1
            advance = 0
            while True:
                num = pileup_bases[base_idx]
                if num.isdigit():
                    advance = advance * 10 + int(num)
                    base_idx += 1
                else:                                                                                                                                                                                                                                                                      
                    break
            base_list.append(base + pileup_bases[base_idx: base_idx + advance])
            base_idx += advance - 1

        elif base == '^':  # start of a read, next character is mapping quality
            base_idx += 1
        # elif base == '$': # end of read with '$' symbol
        base_idx += 1
    base_counter = Counter(base_list)
    return base_counter, _cov

def get_cnt_from_counter(base_counter, _tar_ref, _tar_alt):
    _tar_b = _tar_alt
    if len(_tar_ref) > len(_tar_alt):
        # del
        _n_b = "N" * (len(_tar_ref) - len(_tar_alt))
        _tar_b = "-%s" % (_n_b)
    elif len(_tar_ref) < len(_tar_alt):
        # ins
        _tar_b = "+%s" % (_tar_alt[-(len(_tar_alt) - len(_tar_ref)):])
    return base_counter[_tar_b] if _tar_b in base_counter else 0


def update_af(_bam, _bed, vcf_d_ori, _is_try_update=True, _max_depth=10000, is_log=True):
#     _max_depth = 100000
    _min_mq = 5
    _min_bq = 0
    _samtools_e_flag = 2316
    _bed_region = _bed

    mq_option = ' --min-MQ {}'.format(_min_mq)
    bq_option = ' --min-BQ {}'.format(_min_bq)
    flags_option = ' --excl-flags {}'.format(_samtools_e_flag)
    max_depth_option = ' --max-depth {}'.format(_max_depth)
    bed_option = ' -l {}'.format(_bed) if _bed != "" else ""
    mpileup_cmd = "samtools mpileup {}".format(_bam) + mq_option + bq_option + flags_option + max_depth_option + bed_option
    if is_log:
        print(mpileup_cmd)
    output = subprocess.getoutput(mpileup_cmd)

    all_p = {r[1]: _i for _i, r in enumerate(vcf_d_ori)}
    vcf_d = []
    for row in output.split('\n'):
        row = row.split('\t')
        if len(row) > 3:
            _c, _p, _cov, _s = row[0], int(row[1]), int(row[3]), row[4]
            _s = _s.upper()
            if _p not in all_p:
                continue
            _counter, _cov  = parse_pileup_out(_s)
            _tar_ref, _tar_alt = vcf_d_ori[all_p[_p]][2].upper(), vcf_d_ori[all_p[_p]][3].split(",")[0].upper()
            _n_cnt = get_cnt_from_counter(_counter, _tar_ref, _tar_alt)
            _new_af = 1. * _n_cnt / _cov
                    
            if is_log:
                print("%s %s %s" % (_tar_ref, _tar_alt, _counter), file=sys.stderr)
                print("%s %s %s" % (_n_cnt, _cov, _n_cnt / _cov), file=sys.stderr)
                print("%s" % (vcf_d_ori[all_p[_p]]), file=sys.stderr)
            _new_rec = vcf_d_ori[all_p[_p]][:]
            _new_rec[4] = str(_cov)
            _ori_af = vcf_d_ori[all_p[_p]][-1]
            # vcf indel error, true af is different from called af, try to rescure the indel
            if _is_try_update:
                _new_rec[3] = _tar_alt
                if abs(_new_af - _ori_af) > 0.5 and _new_af < 0.1:
                # try to replace the low af variant with high af
                    print("------\nvariant adjustment", file=sys.stderr)
                    _ori_ref, _ori_alt = _tar_ref, _tar_alt
                    _new_ref, _new_alt = "", ""
                    _new_c = 0
                    try:
                        # only try to replace indel
                        if len(_ori_alt) == len(_ori_ref):
                            raise TypeError()                        
                        _tar_l = [[i, _counter[i]] for i in _counter if i != _ori_ref[0]]
                        _tar_l.sort(key=lambda x: -x[1])
                        _tar_tmp_s, _new_c = _tar_l[0][0], _tar_l[0][1]
                        print("%s, %s %s %s" % (_p, _tar_l, _tar_tmp_s, _new_c), file=sys.stderr)
                
                        if _tar_tmp_s == "*":
                            print("no valid v", file=sys.stderr)
                            _new_ref, _new_alt = "", ""
                            _new_af = 0
                        elif len(_tar_tmp_s) == 1:
                            # snp
                            _new_ref = _ori_ref[0]
                            _new_alt = _tar_tmp_s
                        elif _tar_tmp_s[0] == "+":
                            # ins
                            _new_ref = _ori_ref[0]
                            _new_alt = _new_ref + _tar_tmp_s[1:]
                        else:
                            # del
                            _new_ref = _ori_ref[0] + "N" * (len(_tar_tmp_s) - 1)
                            _new_alt = _ori_ref[0]
                        
                        
                    except:
                        _new_ref, _new_alt = "", ""
#                         pass
                    if _new_ref != "":
                        print("update %s: %s -> %s | %s -> %s" % (_p, _ori_ref, _ori_alt, _new_ref, _new_alt), file=sys.stderr)
                        _new_rec[2] = _new_ref
                        _new_rec[3] = _new_alt
                        _new_af = 1. * _new_c / _cov
                    
            _new_rec.append(_new_af)
            vcf_d.append(_new_rec)
            
            if is_log:
                print("%s\n\n" % (_new_rec), file=sys.stderr)
            
    return vcf_d


def update_vcf(_vcf, _bam, _bed, _out_fn, max_dp_in_check=10000, is_log=True, _is_try_update=True, _is_filter_empty=False):
    print(_vcf, _bam, _bed, _out_fn)
    tree = bed_tree_from(bed_file_path=_bed)
    is_tree_empty = len(tree.keys()) == 0
    vcf_d_ori, v_cnt = read_vcf(_vcf, is_tree_empty, tree, is_snp_only=False)
    
#     _is_try_update = True
    _max_depth=max_dp_in_check
    vcf_d = update_af(_bam, _bed, vcf_d_ori, _is_try_update, _max_depth, is_log)
    
    d = {"%s:%s" % (i[0], i[1]) : i for i in vcf_d}
    with open(_out_fn, 'w') as _O:
        vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (_vcf)))
        for row in vcf_fp.stdout:
            columns = row.strip().split()
            if columns[0][0] == "#":
                _O.write(row)
                continue
            ctg_name, position = columns[0], int(columns[1])
            _id = "%s:%s" % (ctg_name, position)
            if _id not in d:
                if not _is_filter_empty:
                    _O.write(row)
            else:
                if _is_try_update and (columns[3] != d[_id][2] or columns[4] != d[_id][3]):
                    columns[3] = d[_id][2]
                    columns[4] = d[_id][3]
                new_col = "\t".join(columns[:-1])
                if _is_try_update:
                    new_col += "\t%s:%s:%s:%.6f\n" % ("1", ":".join(columns[-1].split(":")[1:-2]), d[_id][4], d[_id][-1])
                else:
                    new_col += "\t%s:%s:%.6f\n" % (":".join(columns[-1].split(":")[:-2]), d[_id][4], d[_id][-1])
                    
                _O.write(new_col)
    
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_fn', type=str, default="input.vcf",
                    help="Truth VCF file input, required")
    
    parser.add_argument('--bam_fn', type=str, default="input.bam",
                        help="input bam file, required")

    parser.add_argument('--bed_fn', type=str, default="input.bed",
                        help="input bed file, optional")
    
    parser.add_argument('--out_fn', type=str, default="output.vcf",
                    help="Truth VCF file input, required")
    
    parser.add_argument('--max_dp_in_check', type=int, default=10000, 
                    help="max coverage scanned for updateing af")
    
    parser.add_argument('--no_try_update', default=False, action='store_true',
                        help="is no try update vcf based on coverage information, optional")
    
    parser.add_argument('--filter_empty', default=False, action='store_true',
                        help="filter no coverage row")
    args = parser.parse_args()
    
    _vcf = args.vcf_fn
    _bam = args.bam_fn
    _bed= args.bed_fn
    _out_fn = args.out_fn
    _max_dp_in_check = args.max_dp_in_check
    _is_try_update = True
    if args.no_try_update:
        _is_try_update = False
    
    _is_filter_empty = False
    if args.filter_empty:
        _is_filter_empty = True
#     print(_is_filter_empty)
#     return 0
#     _is_try_update = False
    
    update_vcf(_vcf, _bam, _bed, _out_fn, _max_dp_in_check, is_log=False, _is_try_update=_is_try_update, _is_filter_empty=_is_filter_empty)

if __name__=='__main__':
    main()
