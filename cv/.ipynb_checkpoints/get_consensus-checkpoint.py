import sys
import os
import re
import shlex
from subprocess import PIPE
from argparse import ArgumentParser
from shared.utils import subprocess_popen, vcf_candidates_from, _run_command
from shared.interval_tree import bed_tree_from, is_region_in
import subprocess
import json
import pandas as pd
from matplotlib import pyplot as plt 
import seaborn as sns
import numpy as np
import hashlib

from cv.update_vcf_af import update_vcf
from cv.parse_hivdb_json import parse_hivdb_json

def read_vcf(vcf_fn):
    vcf_fp = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))
    vcf_list = []
    v_cnt = {'snp': 0, 'indel': 0}
    _header_l = []
    for row in vcf_fp.stdout:
        columns = row.strip().split()
        if columns[0][0] == "#":
            _header_l.append(row.strip())
            continue
        # position in vcf is 1-based
        ctg_name, position = columns[0], columns[1]
        i = columns[:]
        # chr, pos, ref_base, alt_base, qual, info, info, af
#         tar_info = [i[0], int(i[1]), i[3], i[4], i[5], i[-2], i[-1], float(i[-1].split(':')[-1])]
        # chr, pos, ref_base, alt_base, qual, af
        tar_info = [i[0], int(i[1]), i[3], i[4], float(i[-1].split(':')[-1]), columns]
        if len(i[3]) == 1 and all([len(_j) == 1 for _j in i[4].split(',')]) == 1:
            v_cnt['snp'] += 1
        else:
            v_cnt['indel'] += 1
        vcf_list.append(tar_info)
    # sort by chromosome and position
    vcf_list.sort(key=lambda x: (x[0], x[1]))
    print('read %s snp, %s indel from %s' % (v_cnt['snp'], v_cnt['indel'], vcf_fn))
    return vcf_list, _header_l


def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, str):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

_G_map = {'PR': 1799, 'RT': 2096, 'IN': 3776}
def _get_p(_t, _m, _G_map):
    _p = re.findall('(\d+)', _m)[0]
    return _G_map[str(_t)] + 3 * int(_p) - 1

def run_get_consensus(args):
    _all_out_dir = args.out_dir
    _info_tsv = args.tar_tsv
    _ref = args.ref
    _bed = args.bed
    _hivdb_url = args.hivdb_url
    _number_of_read_for_consense = args.number_of_read_for_consense
    _hivdb_url_option = '' if _hivdb_url == '' else "--url %s" % (_hivdb_url)
    _py_s_d = os.path.dirname(os.path.abspath(__file__))
    print(_py_s_d)
    
    # only support one bed region
    with open(_bed, 'r') as F:
        for line in F:
            if len(line.strip()) <= 0:
                continue
            _bed_ctg, _bed_s, bed_e = line.split()
            break
    print(_bed_ctg, _bed_s, bed_e)
    
    cmd = 'mkdir -p %s' % (_all_out_dir)
    _run_command(cmd)
    cmd = 'cp %s %s' % (_info_tsv, _all_out_dir)
    _run_command(cmd)


    _png = '/'.join(_info_tsv.split('/')[:-1]) + '/*.png'
    try:
        cmd = 'cp %s %s' % (_png, _all_out_dir)
        _run_command(cmd)
    except:
        pass
    
    all_c_l = []
    all_c_l.append(['id', 'gene', 'drugClass', 'drugName', 'drugScore', 'resistanceLevel', 'subtype_ori', 'subtype', 'abundance', 'is_in_consensus', 'VAF', 'mutationPos', 'mutation', 'mutationScore', 'mutationType', 'comments'])
    _new_info = []
    _new_info.append(['ID', 'PROPORTION', 'CONSENSUS_FA', 'BAM_CONSENSUS', 'VCF_CONSENSUS', 'VCF_ORI', 'BAM_ORI', ])
    _id = ''
    _id_v_lst = []
    _info_tsv_l = []
    with open(_info_tsv, 'r') as F:
        for line in F:
            if len(line.strip()) <= 0 or line[0] == '_':
                if line[0] == '_':
                    _n_row = line.strip().split()
                    _info_tsv_l.append(_n_row)
                continue
            row = [i for i in line.strip().split() if len(i) > 0]
            _ori_vcf, _ori_bam, _ori_id, _p = row[1], row[2], row[5], row[-2]
            _s_idx = row[0]

            _id = row[4].split('.')[0]
            _s_id = '_'.join(_id.split('_')[1:])
#             _new_id = _ori_id
            _new_id = _s_id + '_' + _s_idx
            _n_row = row
            row[0] = _new_id
            _info_tsv_l.append(_n_row)
            
            _out_dir = "%s/%s" % (_all_out_dir, _new_id)
            cmd = 'mkdir -p %s' % (_out_dir)
            _run_command(cmd)
            
            # vcf from cliar
            _new_vcf = "%s/%s.vcf.gz" % (_out_dir, _new_id)
            _new_bam = "%s/%s.bam" % (_out_dir, _new_id)
            _new_bam_read = "%s/%s_ori_r.fasta" % (_out_dir, _new_id)
            
            _new_cs_dir = "%s/%s_flye" % (_out_dir, _new_id)
            # consensus fasta, bam, vcf
            _new_cs = "%s/%s.fasta" % (_all_out_dir, _new_id)
            _new_cs_bam = "%s/%s_cs.bam" % (_out_dir, _new_id) 
            _new_cs_vcf_tmp = "%s/%s_cs_tmp.vcf" % (_out_dir, _new_id) 
            _new_cs_vcf = "%s/%s_cs.vcf" % (_out_dir, _new_id) 
            print(_ori_vcf, _ori_bam, _ori_id, _p, _id, _s_id, _s_idx)
            
            # get bam
            cmd = 'cp %s %s; cp %s.bai %s.bai;' % (_ori_bam, _new_bam, _ori_bam, _new_bam)
            _run_command(cmd)
            os.system(cmd)
            
            # get vcf
            cmd = 'bcftools filter -e "FORMAT/AF<=0.2" %s | bcftools view -O z -o %s; bcftools index %s' % (_ori_vcf, _new_vcf, _new_vcf)
            _run_command(cmd)
            
#             get consensus
#             get read from bam
            cmd = 'samtools fasta %s > %s_1' % (_new_bam, _new_bam_read)
            _run_command(cmd)
            cmd = 'head -n %d %s_1 > %s; rm %s_1' % (_number_of_read_for_consense, _new_bam_read, _new_bam_read, _new_bam_read)
            _run_command(cmd)
            
            # run flye
            cmd = 'flye --nano-raw %s --threads 32 --out-dir %s -m 1000 -g 5k' % (_new_bam_read, _new_cs_dir)
            _run_command(cmd)

            # run alignmnet for visulization
            cmd = 'minimap2 -ax map-ont %s %s/assembly.fasta | samtools sort | samtools view -F 2048 -b > %s; samtools index %s' % (_ref, _new_cs_dir, _new_cs_bam, _new_cs_bam)
            _run_command(cmd)
            
            # generate consensus for visulization
            cmd = "samtools mpileup -f %s %s | python %s/mpileup2vcf.py --sample_n %s | bcftools view > %s" % (_ref, _new_cs_bam, _py_s_d, _new_id, _new_cs_vcf_tmp)
            _run_command(cmd)

            # get consensus varaints and AF
            update_vcf(_new_cs_vcf_tmp, _new_bam, _bed, _new_cs_vcf, max_dp_in_check=10000, is_log=False)
            
            _cs_vcf_lst, _ = read_vcf(_new_cs_vcf)
            print(_cs_vcf_lst[0])
            # flatten varaints list
            _tar_v_list = ','.join(["%s:%s,%s" % (_tmp_v[1], _tmp_v[2], _tmp_v[3]) for _tmp_v in _cs_vcf_lst])
            _hash_m = hashlib.sha256(_tar_v_list.encode('utf-8')).hexdigest()
            print(_new_id, _p, _hash_m)
            _is_same_cns = False
            _new_id_v_lst = []
            for _rec in _id_v_lst:
                _tmp_s, _tmp_p, _tmp_h = _rec[0], _rec[1], _rec[2]
                if _tmp_h == _hash_m:
                    print('>>>warning, same consensus found in %s, %s match %s' % (_id, _tmp_s, _new_id), file=sys.stderr)
                    _tmp_p = float(_tmp_p) + float(_p)
                    _is_same_cns = True
                _new_id_v_lst.append([_tmp_s, _tmp_p, _tmp_h])
            if not _is_same_cns:
                _new_id_v_lst.append([_new_id, _p, _hash_m])
            _id_v_lst = _new_id_v_lst[:]
            if _is_same_cns:
                continue
                
            cmd = 'echo ">%s" > %s; tail -n+2 %s/assembly.fasta >> %s' % (_new_id, _new_cs, _new_cs_dir, _new_cs)
            _run_command(cmd)

            # get drug resistance report
            cmd = 'sierrapy fasta %s %s > %s/%s_report.json' % (_new_cs, _hivdb_url_option, _all_out_dir, _new_id)
            _run_command(cmd)

            _f_p = '%s/%s_report.json' % (_all_out_dir, _new_id)
            _f_o = '%s/%s_report.tsv' % (_all_out_dir, _new_id)
            
            all_l = []
            all_l.append(['gene', 'drugClass', 'drugName', 'drugScore', 'resistanceLevel', 'subtype_ori', 'subtype', 'abundance', 'is_in_consensus', 'VAF', 'mutationPos', 'mutation', 'mutationScore', 'mutationType', 'comments'])
            
            def get_v_from_vcf_in_range(_vcf_lst, _l_p, _r_p):
                tar_v_l = []
                for i in range(len(_vcf_lst)):
                    _v_rec = _vcf_lst[i]
                    _v_p = _v_rec[1]
                    if _v_p < _l_p:
                        continue
                    if _v_p > _r_p:
                        break
                    tar_v_l.append(_v_rec[:5]) 
                return tar_v_l

            with open(_f_p, 'r') as f:
                _data = json.load(f)
                tmp_all_l = parse_hivdb_json(_data)
                if len(tmp_all_l) > 1:
                    for i in tmp_all_l[1:]:
                        _gene_n, _mutation = i[0], i[5]
                        _type = i[4]
#                         if _type == "Susceptible":
#                             continue
                        _pos = _get_p(_gene_n, _mutation, _G_map)
#                         all_l.append([i['gene']['name'],  all_mut, j['drugClass']['name'], j['drug']['name'], j['text'],  _ori_id, _new_id, _p, '1', _pos,])
                        
                        _tmp_tar_v_l = get_v_from_vcf_in_range(_cs_vcf_lst, _pos-2, _pos)
                        _tar_v_lst_s = "|".join(["%s,%s,%s,%s" % (_tv[1], _tv[2], _tv[3], _tv[4]) for _tv in _tmp_tar_v_l])
                        _tar_v_lst_p = 1 if len(_tmp_tar_v_l) <= 0 else min([_tv[4] for _tv in _tmp_tar_v_l]) 
#                         print(_tar_v_lst_s, _tar_v_lst_p)
                        all_l.append(i[:5] + [_ori_id, _new_id, '%.5f' % float(_p), '1', '%s|%s' % (_tar_v_lst_p, _tar_v_lst_s), _pos] + i[5:])
                        all_c_l.append([_id] + all_l[-1])
            
            # check low af report  
            _py_s_d = os.path.dirname(os.path.abspath(__file__))
            print('checking low af variant')
            cmd = "python %s/../cv.py get_low_v_db --cns_vcf %s --ori_vcf %s --sample_ori_id %s --sample_id %s --vcf_o_dir %s --ref %s" % \
                (_py_s_d, _new_cs_vcf, _new_vcf, _ori_id, _new_id, _out_dir, _ref)
            _run_command(cmd)
                
            _low_af_rst_f = "%s/%s_report_low_af.tsv" % (_out_dir, _new_id)
            if os.path.exists(_low_af_rst_f):
                # read low af relusts
                with open(_low_af_rst_f, 'r') as F:
                    for line in F:
                        line = line.strip()
                        if len(line) <= 0 or line[0] == 'g':
                            continue
                        row = line.split('\t')
                        row[7] = _p
                        print('hit low af', row)
#                         row = row[:6] + [_ori_id] + row[6:]
                        all_l.append(row)
                        all_c_l.append([_id] + all_l[-1])
            print('done check low af')
            
            with open(_f_o, 'w') as F:
                for l in all_l:
                    F.write("%s\n" % ("\t".join([str(k) for k in l])))

            
            _new_info.append([_new_id, _p, _new_cs, _new_cs_bam, _new_cs_vcf, _new_vcf, _new_bam])
    
    # update proportion with corrected info.
    def _update_p(new_p_l, tar_l, _id_idx, _p_idx):
        # id: proportion mapping
        _m = {i[0]: i[1] for i in new_p_l}
        _new_l = []
        for _i in range(len(tar_l)):
            r = tar_l[_i]
            if _i == 0:
                _new_l.append(r)
                continue
            if r[_id_idx] not in _m:
                continue
            r[_p_idx] = _m[r[_id_idx]]
            _new_l.append(r)
        return _new_l
    print(_id_v_lst)
    _u_info_tsv_l  = _update_p(_id_v_lst, tar_l = _info_tsv_l, _id_idx = 0, _p_idx = -2)
    _info_tsv_l = _u_info_tsv_l
    _u_new_info  = _update_p(_id_v_lst, tar_l = _new_info, _id_idx = 0, _p_idx = 1)
    _new_info = _u_new_info
    _u_all_c_l  = _update_p(_id_v_lst, tar_l = all_c_l,  _id_idx = 7, _p_idx = 8)
    all_c_l = _u_all_c_l
    
    # subtype's info.
    new_tsv = "%s/all_info.tsv" % (_all_out_dir)
    with open(new_tsv, 'w+') as F:
        for r in _info_tsv_l:
            F.write("\t".join([str(k) for k in r]))
            F.write("\n")
            
    # consensus info. 
    new_tsv = "%s/info.tsv" % (_all_out_dir)
    with open(new_tsv, 'w+') as F:
        for r in _new_info:
            F.write("\t".join([str(k) for k in r]))
            F.write("\n")
            
    # drug resistance report
    new_tsv = "%s/all_report.tsv" % (_all_out_dir)
    with open(new_tsv, 'w+') as F:
        for l in all_c_l:
            F.write("%s\n" % ("\t".join([str(k) for k in l])))
#             F.write("%s\n" % ("\t".join(l)))
    print('all report at %s' % (new_tsv))
    v_all_c_l = [i for i in all_c_l if i[2] != 'NA']
    if len(v_all_c_l) <= 1:
        print('no drug resistance found')
        return 0
#     return
    _sample_id = _id
    _tar_f = "%s/all_report.tsv" % (_all_out_dir)

    _all_df = pd.read_csv(_tar_f, sep='\t')
#     _all_df = _all_df[_all_df['drugName'] != 'NA']
    _all_df = _all_df.dropna()
    if len(_all_df) < 1:
        print('no drug resistance found')
        return 0
    new_p_map = {'High-Level Resistance': 'HL',
                 'Intermediate Resistance': 'IL',
                 'Low-Level Resistance': 'LL',
                 'Potential Low-Level Resistance': 'PLL',
                 'Susceptible': 'S'}
    _all_df['resistance_level'] = _all_df.apply(lambda x: new_p_map[x['resistanceLevel']], axis = 1)
    _all_df['subtype'] = _all_df.apply(lambda x: x['subtype'].split('_')[-1], axis = 1)
    _all_df = _all_df[_all_df['resistance_level'] != "S"]

    if len(_all_df) > 0:
    #     fig = plt.figure(figsize = (12, 6))
        fig, (ax2, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 3]}, figsize = (10, 5))
        fig.suptitle('Drug resistance mutation in %s' % (_sample_id))
        _tmp_df = _all_df[['gene', 'mutation', 'abundance', 'subtype']]
        _tmp_df = _tmp_df.drop_duplicates()
        tar_dr = _tmp_df.groupby(['gene', 'mutation', 'subtype'])['abundance'].sum().unstack(["subtype"]).fillna(0)
        tar_dr.plot.bar(stacked=True, ax = ax2)

        ax2.set_ylabel("subtype's abudance")
        ax2.set_ylim([0, 1.05])
        ax2.xaxis.grid(False, which='both')

        ax2.margins(.05)
        for label in ax2.get_xticklabels():
            label.set_rotation(90)
        ax2.get_legend().remove()

        _tmp_df = _all_df[['resistance_level', 'drugName', 'abundance', 'subtype']]
        _tmp_df = _tmp_df.drop_duplicates()
        tar_dr = _tmp_df.groupby(['resistance_level', 'drugName', 'subtype'])['abundance'].sum().unstack(["subtype"]).fillna(0)
        tar_dr.plot.bar(stacked=True, ax=ax1)

        ax1.set_ylim([0, 1.05])

        ax1.margins(.05)
        for label in ax1.get_xticklabels():
            label.set_rotation(90)
        ax1.legend(loc='center left', bbox_to_anchor=(1.02, 0.6), title="subtype",fontsize='x-small', fancybox=True)
        ax1.xaxis.grid(False, which='both')

    #     plt.tight_layout()
        plt.savefig('%s/all_mutation_and_drug_ressitant_barplot.png' % (_all_out_dir), dpi=100, bbox_inches = "tight")
        
        _tmp_df = _all_df[['gene', 'mutation', 'abundance', 'subtype']]
        _tmp_df = _tmp_df.drop_duplicates()
        tar_dr = _tmp_df.groupby(['gene', 'mutation', 'subtype'])['abundance'].sum().unstack(["subtype"]).fillna(0)
        _tar = tar_dr.transpose()
        _tar.loc['All'] = _tar.sum(0)
        _d, _s = _tar.shape
        _tar = _tar.replace(0, np.nan)

        _tmp_df = _all_df[['resistance_level', 'drugName', 'abundance', 'subtype']]
        _tmp_df = _tmp_df.drop_duplicates()
        tar_dr1 = _tmp_df.groupby(['resistance_level', 'drugName', 'subtype'])['abundance'].sum().unstack(["subtype"]).fillna(0)

        _tar1 = tar_dr1.transpose()
        _tar1.loc['All'] = _tar1.sum(0)

        _d1, _s1 = _tar1.shape
        print(_d1, _s1)
        _tar1 = _tar1.replace(0, np.nan)

        _cell_w = 1
        fig, axes = plt.subplots(1, 2, figsize=((_s + _s1 + 3.) * _cell_w, (_d + 1) * _cell_w), gridspec_kw={'width_ratios': [(_s) * _cell_w, (_s1+ 2) * _cell_w]})

        ax1 = axes[0]
        ax1 = sns.heatmap(
            _tar, 
            vmin=0, vmax=1, center=0.7,
            cmap=sns.light_palette("firebrick", as_cmap=True),annot=True, fmt='.2f',
            square=True, cbar_kws={'label': 'abundance', "shrink": .87}, linewidths=0.2, ax=ax1, cbar=False
        )
        plt.title("Linear graph")

        ax1.set_xticklabels(
            ax1.get_xticklabels(),
            rotation=45,
            horizontalalignment='center'
        );
        ax1.set_title('')
        ax1.tick_params(left=False)

        for _, spine in ax1.spines.items():
            spine.set_visible(True)

        # fig, ax2 = plt.subplots(figsize=(_s * 0.9 + 0.5, _d * 0.9 + 0.5))
        ax2 = axes[1]
        ax2 = sns.heatmap(
            _tar1, 
            vmin=0, vmax=1, center=0.7,
            cmap=sns.light_palette("firebrick", as_cmap=True),annot=True, fmt='.2f',
            square=True, cbar_kws={'label': 'abundance', "shrink": .7}, linewidths=0.2, ax=ax2, cbar=True
        )

        ax2.set_xticklabels(
            ax2.get_xticklabels(),
            rotation=45,
            horizontalalignment='center'
        );    
        ax2.tick_params(left=False)
        ax2.set(ylabel=None)

        for _, spine in ax2.spines.items():
            spine.set_visible(True)
        ax1.annotate("*Resistance level: HL: High-Level Resistance, IL: Intermediate Resistance,\n                              LL: Low-Level Resistance, PLL: Potential Low-Level Resistance", 
                    size=11, xy = (-0.1, -2.5 / (_d * 0.9 + 0.5)), xycoords='axes fraction')
        ax2.set_title('Heatmap of gene and drug resistance(s) in subtypes                            ', pad=20)
        # ax2.set_title('')
        # fig.suptitle('Heatmap of gene and drug resistance(s) in subtypes')
        plt.savefig('%s/all_mutation_and_drug_ressitant_heatmap.png' % (_all_out_dir), dpi=100, bbox_inches = "tight")
        


    return 0
    
def main():
    parser = ArgumentParser(description="gather results files and get consensus")

    parser.add_argument('--out_dir', type=str, default="output_dir",
                    help="Output folder, required")

    parser.add_argument('--tar_tsv', type=str, default="tar.tsv",
                        help="tsv file from clustering results")
    
    parser.add_argument('--ref', type=str, default="ref.fa",
                        help="ref fasta file")
    
    parser.add_argument('--bed', type=str, default="region.bed",
                        help="bed region")
    
    parser.add_argument('--number_of_read_for_consense', type=int, default=1000,
                        help="number of original read for generating consense")
    
    parser.add_argument('--hivdb_url', type=str, default="",
                        help="hivdb url defalut query from internet, for localize the HIVDB, please check https://github.com/hivdb/sierra, and update this setting accordingly, e.g. \
                        by using --hivdb_url http://localhost:8111/sierra/rest/graphql")
    
    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    run_get_consensus(args)

if __name__ == "__main__":
    main()
