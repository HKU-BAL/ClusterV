import argparse
from argparse import ArgumentParser
import sys
import os
import re
import shlex
from shared.utils import subprocess_popen, _run_command
import math
import json
from cv.parse_hivdb_json import parse_hivdb_json


def read_ref(ref_fn):
    ref_d = {}
    _ct, _s = '', ''
    with open(ref_fn, 'r') as F:
        for line in F:
            line = line.strip()
            if len(line) <= 0:
                continue
            if line[0] == '>':
                if _ct != '':
                    ref_d[_ct] = _s
                _ct = line[1:].split(' ')[0]
                _s = ''
            else:
                _s += line
    if _ct != '':
        ref_d[_ct] = _s
    return ref_d


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
    
# check if have partial overalp
def _check_olp(p1, p2):
    _olp_1 = (p1[0] >= p2[0] and p1[0] <= p2[1]) and not (p1[1] >= p2[0] and p1[1] <= p2[1])
    _olp_2 = not (p1[0] >= p2[0] and p1[0] <= p2[1]) and (p1[1] >= p2[0] and p1[1] <= p2[1])
    return _olp_1 or _olp_2

# get varaints only at _v1 but not at _v2, also filter ill defined
def get_isec_10(_v1, _v2):
    _idx1, _idx2 = 0, 0
    _tar_l = []
    while _idx1 < len(_v1):
        _cur_r_pos = _v1[_idx1][1]
        _cur_ref_l = len(_v1[_idx1][2])
        _cur_cov = (_cur_r_pos, _cur_r_pos + _cur_ref_l - 1)
        while _idx2 < len(_v2) and _v2[_idx2][1] < _cur_cov[0]:
            _idx2 += 1
        if _idx2 < len(_v2) and _check_olp(_cur_cov, (_v2[_idx2][1], _v2[_idx2][1] + len(_v2[_idx2][2]) - 1)):
            # have olp
#             print('match %s %s'  % (_v1[_idx1], _v2[_idx2]))
            pass
        elif _idx2 >=  len(_v2) or _v2[_idx2][1] > _cur_cov[1]:
#             print('hit %s %s' % (_v1[_idx1], _v2[_idx2]))
            _tar_l.append(_v1[_idx1])
        _idx1 += 1
    return _tar_l
  

codon_map = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

protein_s = {'PR': 'PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF', 
             'RT': 'PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDKDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYAGIKVKQLCKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLKTGKYARMRGAHTNDVKQLTEAVQKIATESIVIWGKTPKFKLPIQKETWEAWWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRETKLGKAGYVTDRGRQKVVSLTDTTNQKTELQAIHLALQDSGLEVNIVTDSQYALGIIQAQPDKSESELVSQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLVSAGIRKVL',
             'IN': 'FLDGIDKAQEEHEKYHSNWRAMASDFNLPPVVAKEIVASCDKCQLKGEAMHGQVDCSPGIWQLDCTHLEGKIILVAVHVASGYIEAEVIPAETGQETAYFLLKLAGRWPVKTIHTDNGSNFTSTTVKAACWWAGIKQEFGIPYNPQSQGVVESMNKELKKIIGQVRDQAEHLKTAVQMAVFIHNFKRKGGIGGYSAGERIVDIIATDIQTKELQKQITKIQNFRVYYRDSRDPLWKGPAKLLWKGEGAVVIQDNSDIKVVPRRKAKIIRDYGKQMAGDDCVASRQDED'}

protein_l = {'PR': (1799, 1799 + 3 * len(protein_s['PR'])), 'RT': (2096, 2096 + 3 * len(protein_s['RT'])), 'IN': (3776, 3776 + 3 * len(protein_s['IN']))}

def get_codon(_c, codon_map):
    _s = ''
    for i in _c:
        _s = _s + (i if i != 'T' else 'U')
    _codon = codon_map[_s]
    _codon = _codon if _codon != "STOP" else "*"
    return _codon


# dna mutation to amino acid mutation
def bp_to_aa(m, ref_d, codon_map, protein_l):
    m_ctg, m_pos, m_ref, m_alt = m[0], int(m[1]), m[2], m[3].split(',')[0]
    for g in protein_l:
        _g_r = protein_l[g]
        if m_pos > _g_r[0] - 1 and m_pos < _g_r[1] - 1:
            # hit
            _codon_s = math.floor((m_pos - _g_r[0]) / 3) * 3 + _g_r[0]
            _protein_p = math.floor((m_pos - _g_r[0]) / 3) + 1
#             _codo_s = math.floor((m_pos) / 3) * 3 + 1
            _ori_codon = ref_d[m_ctg][_codon_s - 1: _codon_s - 1 + 3]
            _ori_codon_offset = m_pos - _codon_s
            _ref_code = get_codon(_ori_codon, codon_map)
            print(m_pos, m_ref, m_alt, g, _protein_p, _ori_codon, _ori_codon_offset, _ref_code, end=', ')
            print('True %s' % protein_s[g][_protein_p - 1], end=', ')
            if len(m_ref) == len(m_alt):
                # snp
                _new_codon = ''.join([c if i != _ori_codon_offset else m_alt for i, c in enumerate(_ori_codon)])
                _new_c = get_codon(_new_codon, codon_map)
            elif len(m_ref) > len(m_alt):
                # del
                _new_codon = ''
                _new_c = 'del'
            else:
                _new_codon = ''
                _new_c = 'ins'
            print(_new_codon, _new_c, end=', ')
            _n_ref_code = protein_s[g][_protein_p - 1]
            _tar_aa = "%s:%s%s%s" % (g, _n_ref_code, str(_protein_p), _new_c)
            print(_tar_aa)
            if _n_ref_code != _new_c:
                return _tar_aa
        else:
            continue

    return ''

def write_to_vcf(_out_fn, _header, _v):
    with open(_out_fn, 'w') as _O:
        for i in _header:
            _O.write("%s\n" % (i))
        for i in _v:
            _O.write("%s\n" % ("\t".join(i[-1])))

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

def run(args):
    _ori_vcf = args.ori_vcf # all v vcf from clair-ensemble
    _cns_vcf = args.cns_vcf # consensus vcf from flye
    _sample_id = args.sample_id # sample_id
    _ori_sample_id = args.sample_ori_id # ori_sample_id
    _vcf_o_dir = args.vcf_o_dir # consensus records dir from flye
    _low_af_vcf = "%s/%s_low_af.vcf" % (_vcf_o_dir, _sample_id)
    _hivdb_url = args.hivdb_url
    _hivdb_url_option = '' if _hivdb_url == '' else "--url %s" % (_hivdb_url)
    _ref = args.ref

#     _cns_vcf = "%s/%s_cns.vcf" % (_vcf_o_dir, _sample_id)
#     _py_s_d = os.path.dirname(os.path.abspath(__file__))
#     print(_py_s_d)
# #     # get consensus vcf
#     cmd = "samtools mpileup -f %s %s | python %s/mpileup2vcf.py | bcftools view > %s" % (_ref, _cns_bam, _py_s_d, _cns_vcf)
#     _run_command(cmd)

    _ori_v, _header = read_vcf(_ori_vcf)
    _cns_v, _ = read_vcf(_cns_vcf)
    _tar_v_ori = get_isec_10(_ori_v, _cns_v)
    write_to_vcf(_low_af_vcf, _header, _tar_v_ori)
    _tar_v = [i[:-1] for i in _tar_v_ori]
    ref_d = read_ref(_ref)
#     print([(i, len(ref_d[i])) for i in ref_d])
    _f_log = "%s/%s_low_af_lst.log" % (_vcf_o_dir, _sample_id)
    print('log file at %s' % (_f_log))
    _F = open(_f_log, 'w')

    _all_aa = {}
    for _m in _tar_v:
        _aa = bp_to_aa(_m, ref_d, codon_map, protein_l)
        _F.write("%s %s\n" % (_m, _aa))
        if _aa == '':
            continue
#         _all_aa.append(_m +  [_aa])
        _all_aa[_aa] = _m
    _F.close()
    if len(_all_aa) <= 0:
        print('no resistance mutation in low af')
        return 0
    print(_all_aa)
    print(' '.join([i for i in _all_aa]))
    _tar_seq = ' '.join([i for i in _all_aa])
    cmd = "sierrapy mutations %s %s > %s/%s_report_low_af.json" % (_tar_seq, _hivdb_url_option, _vcf_o_dir, _sample_id)
    _run_command(cmd)

    # parse report
    _f_p = "%s/%s_report_low_af.json" % (_vcf_o_dir, _sample_id)
    _f_o = "%s/%s_report_low_af.tsv" % (_vcf_o_dir, _sample_id)
    all_l = []
    all_l.append(['gene', 'drugClass', 'drugName', 'drugScore', 'resistanceLevel', 'subtype_ori', 'subtype', 'abundance', 'is_in_consensus', 'VAF', 'mutationPos', 'mutation', 'mutationScore', 'mutationType', 'comments'])
    with open(_f_p, 'r') as f:
        _data = json.load(f)
        tmp_all_l = parse_hivdb_json(_data)
        if len(tmp_all_l) > 1:
            for i in tmp_all_l[1:]:
                _gene_n, _mutation = i[0], i[5]
                _type = i[4]
#                 if _type == "Susceptible":
#                     continue
                _pos = _get_p(_gene_n, _mutation, _G_map)
                _tmp_all_mut = "%s:%s" % (_gene_n,  _mutation)
                _low_af_p = 0 if _tmp_all_mut not in _all_aa else _all_aa[_tmp_all_mut][-1]
                all_l.append(i[:5] + [_ori_sample_id, _sample_id, '0', '0', str(_low_af_p), _pos] + i[5:])

    all_v_l = [i for i in all_l if i[1] != 'NA']
    if len(all_v_l) > 1:
        print('drug resitance found in low af -------------------')
        for l in all_v_l:
            print(l)

    with open(_f_o, 'w') as F:
        for l in all_l:
            F.write("%s\n" % ("\t".join([str(k) for k in l])))



def main():
    parser = ArgumentParser(description="get low AF variants and check HIVDB db")

    parser.add_argument('--ori_vcf', type=str, default='ori_vcf.vcf')
    parser.add_argument('--cns_vcf', type=str, default='cns.vcf')
    parser.add_argument('--sample_id', type=str, default='sample_id')
    parser.add_argument('--sample_ori_id', type=str, default='ori_id')
    parser.add_argument('--vcf_o_dir', type=str, default='vcf_o_dir')
    parser.add_argument('--ref', type=str, default='ref.fasta')
    parser.add_argument('--hivdb_url', type=str, default="",
                        help="hivdb url defalut query from internet")
    args = parser.parse_args()

    run(args)

if __name__=='__main__':
    main()
