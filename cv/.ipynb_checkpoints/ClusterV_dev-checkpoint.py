import sys
import os
import shlex
from subprocess import PIPE
from argparse import ArgumentParser
from shared.utils import subprocess_popen, vcf_candidates_from
from shared.interval_tree import bed_tree_from, is_region_in
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import pysam
from scipy.cluster.hierarchy import ward, fcluster
from scipy.signal import find_peaks
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
import subprocess
import random

def _reduce_backgroup(h):
    _background_dis = np.array([0.00000000e+00, 2.18483723e-04, 1.96635351e-03, 6.55451169e-04,
       6.55451169e-04, 6.55451169e-04, 6.55451169e-04, 1.74786978e-03,
       1.74786978e-03, 8.73934892e-04, 4.36967446e-04, 2.18483723e-03,
       1.96635351e-03, 2.62180468e-03, 3.27725584e-03, 4.58815818e-03,
       5.89906052e-03, 1.00502513e-02, 2.88398514e-02, 9.59143544e-02,
       8.35044789e-01, 0.00000000e+00])
    _background_dis_w = (1-_background_dis)**2
    _norm_h = h * _background_dis_w
    return _norm_h

# find AF peaks from AF list
def find_peaks_af(var_list, out_dir_n='.', _sample_id='peaks', _sample_idx='1', _end_bnd=0.95, is_check_snp_only=True):
    # number of bins in histgram
    n_bins = 20
    peaks_height = 0.1
#     peaks_height = 0.03
    peaks_height = 0.03
#     _end_bnd = 0.95

    af_arr = []
    if is_check_snp_only:
        # get all snp af
        af_arr = [float(i[-1]) for i in var_list if len(i[2])==1 and all([len(_k) == 1 for _k in i[3].split(',')])]
    else:
        af_arr = [float(i[-1]) for i in var_list]
    af_arr = [i for i in af_arr if i < _end_bnd]
    if len(af_arr) < 1:
        return [], 0
#     h = plt.hist(af_arr, bins=n_bins, density=True);
    n_bins_arr = np.arange(0, n_bins+1) / n_bins
#     print(plt.hist(af_arr, bins=n_bins_arr))
    h = plt.hist(af_arr, bins=n_bins_arr, density=True);
    plt.clf()

    density = h[0]
    values = h[1][:-1] + np.diff(h[1])[0] / 2
    
    # add padding to get peak in borders
    density = np.insert(density, 0, 0)
    density = np.insert(density, len(density), 0)
    values = np.insert(values, 0, 0 - 0.05)
    values = np.insert(values, len(values), 1 + 0.05)

    fig = plt.figure(figsize = (6, 5))
    ax = fig.add_subplot(1,1,1)
    norm_density = density / density.sum()
#     print(norm_density)
    norm_density = _reduce_backgroup(norm_density)
    norm_density = norm_density / norm_density.sum()
#     print(norm_density)
    weights = np.ones_like(af_arr) / (len(af_arr))
#     plt.hist(af_arr, bins=n_bins,  alpha=.8, weights=weights)
    plt.hist(af_arr, bins=n_bins_arr,  alpha=.8, weights=weights)
    
    required_dis = max(2, int(np.ceil(0.10 / (1. / n_bins)) + 1))
    threshold = 0.005
    peaks = find_peaks(norm_density, height=peaks_height, threshold=threshold, distance=required_dis)[0]
#     peaks = find_peaks(norm_density, height=peaks_height, distance=required_dis)[0]
    if len(peaks) < 1:
        return [], 0
    plt.plot(values, norm_density);
    for peak in peaks:
        plt.axvline(values[peak], color='r')
        plt.text(values[peak], 0.01, ' %.2f' % (values[peak]))
    plt.xlim([0, _end_bnd+0.05])

    plt.xlabel('AF')
    plt.ylabel('density')
    ax.title.set_text('%s %s %s, variants\' AF histogram plot' % (_sample_id, _sample_idx, _end_bnd))
    plt.savefig('%s/hist_plot_%s.png' % (out_dir_n, _end_bnd))
#     plt.clf()
    rst_peaks = []
    for _i in range(len(peaks)):
#         _tmp_af = min(0.92, values[peaks][_i])
        _tmp_af = values[peaks][_i]
        if _tmp_af < 0.95:
            rst_peaks.append([_tmp_af, norm_density[peaks][_i]])
    if len(rst_peaks) < 1:
        return [], 0
    rst_peaks.sort(key=lambda x: -x[1])
    max_peak_af = rst_peaks[0][0]
    print('number of peaks: %s, maximum AF for peak is: %.2f' % (len(peaks), max_peak_af))

    return rst_peaks, max_peak_af

    
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
        # print(ctg_name, position)
        if not (is_tree_empty or is_region_in(tree, ctg_name, int(position))):
            #print('out of bed pos %s' % (position))
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

def update_af(_bam, _bed, vcf_d):
    _max_depth = 100000
    mpileup_cmd = 'samtools mpileup -q 5 -Q 0 -d %d %s -l %s' % (_max_depth, _bam, _bed)
    if _bed == None:
        mpileup_cmd = 'samtools mpileup -q 5 -Q 0 -d %d %s' % (_max_depth, _bam)
    output = subprocess.getoutput(mpileup_cmd)
    # print(output)

    all_l = []
    def count_s(s):
        D = {'A':0, 'C':0, 'G': 0, 'T':0}
        for i in s:
            if i in D:
                D[i] += 1
        return D

    all_p = {r[1]: _i for _i, r in enumerate(vcf_d)}
    for row in output.split('\n'):
        row = row.split('\t')
        if len(row) > 3:
            _c, _p, _cov, _s = row[0], int(row[1]), int(row[3]), row[4]
            if _p not in all_p:
                continue
            _tar_c = count_s(_s.upper())[vcf_d[all_p[_p]][3]]
            _af = _tar_c * 1.0 / _cov
            vcf_d[all_p[_p]].append(_af)
    return vcf_d

# sort variants by quality and return top N variants for clustering feature
# select snp only
# return the top N variants
def select_top_v(vcf_list, MAX_SNP_IN_CLUSTER=10):
#     MAX_SNP_IN_CLUSTER = 10
    tmp_vcf_d = vcf_list[:]
    tmp_vcf_d.sort(key=lambda x: (-int(x[4])))
    tmp_vcf_d = [i for i in tmp_vcf_d if len(i[2]) == 1 and len(i[3]) == 1]
    tar_vcf_d = tmp_vcf_d[:MAX_SNP_IN_CLUSTER]
    return tar_vcf_d

def basequal_to_prob(bq):
    if bq=='':
        return bq
    else:
        return 1-pow(10,bq/-10)

# # # encode base with basequality
# def transform_nucl_bq(x,bq):
#     print(bq)
#     if x=='A':
#         return [bq,0,0,0]
#     elif x=='T':
#         return [0,bq,0,0]
#     elif x=='C':
#         return [0,0,bq,0]
#     elif x=='G':
#         return [0,0,0,bq]
#     elif x=='DEL':
#         return [0,0,0,0]
    
# encode base with basequality
def transform_nucl_bq(x,bq):
#     bq = 0 if x=='DEL' else (5 if bq < 5 else 30)
    bq = 30
    new_l = [0] * 5
    if x=='A':
        new_l[0] = bq
    elif x=='T':
        new_l[1] = bq
    elif x=='C':
        new_l[2] = bq
    elif x=='G':
        new_l[3] = bq
    elif x=='DEL':
        new_l[4] = bq
    return new_l
    
# encode base with basequality
def transform_nucl_bq_r(x,bq,ref):
    if x==ref:
        return [bq]
    else:
        return [0]

# get read base and quality information from [alignment and selected variants sites]
# return a ditionary with read_id and encoded base
def get_read_v_matrix(bam_fn, tar_vcf_d, MAX_DEPTH=150000):
#     print('get read count matrix.')
    readid_dict=dict()
    read_cov_dict = []
    with pysam.AlignmentFile(bam_fn, "rb" ) as samfile:
        for entry in tar_vcf_d:
            pos=int(entry[1])
            ctg_name = entry[0]
            ref=entry[2]
            alt=entry[3]
#             print(ctg_name, pos, ref, alt)
            _coverage = 0
#             tmp_bq = []
            atcg_cnt = {'A': 0, 'C': 0, 'G': 0, 'T':0, 'DEL': 0}
            for pileupcolumn in samfile.pileup(ctg_name,pos-1,pos,truncate=True,min_base_quality=0,min_mapping_quality=5,flag_filter=2316,max_depth=MAX_DEPTH):
                _coverage = pileupcolumn.n
                for pileupread in pileupcolumn.pileups:
                    bq = 0
                    if pileupread.is_del or pileupread.is_refskip:
                        base='DEL'
                    else:
                        base=pileupread.alignment.query_sequence[pileupread.query_position]
                        bq=pileupread.alignment.query_qualities[pileupread.query_position]
                    if pileupread.alignment.query_name not in readid_dict:
                        readid_dict[pileupread.alignment.query_name]=dict()
#                     tmp_bq.append(bq)
                    readid_dict[pileupread.alignment.query_name][f'{pileupcolumn.pos+1} DP:{_coverage} {ref}>{alt}']=transform_nucl_bq(base,bq)
#                     readid_dict[pileupread.alignment.query_name][f'{pileupcolumn.pos+1} DP:{pileupcolumn.n} {ref}>{alt}']=(base,bq)
                    atcg_cnt[base] += 1
#             _af = 1.0 * atcg_cnt[alt] / _coverage
            _af = 0
            read_cov_dict.append([ctg_name, pos, ref, alt, _coverage, _af, atcg_cnt])
    return readid_dict, read_cov_dict
        
def _get_sample_times(d, s_t=5, MAX_SNP_IN_CLUSTER=25, _step=3):
    _down_l = []
    _tmp_c = min(len(d), MAX_SNP_IN_CLUSTER )
    if len(d) > MAX_SNP_IN_CLUSTER + 10:
        print('ramdom sampling from candidates')
        # random sample
        _run_i = 0
        while _run_i < s_t:
            yield d[:_tmp_c]
            random.shuffle(d)
            _run_i += 1
    else:
        print('down sampling from candidates')
        # downsample
        _run_i = 0
        s_t = min(s_t, int(len(d) / 5 * 1.5))
        while _run_i < s_t and _tmp_c >= 5:
            yield d[:_tmp_c]
            _run_i += 1
            _tmp_c -= _step
            if _tmp_c < 5:
                random.shuffle(d)
                _tmp_c = 5
    yield []
    

def get_clusers_from_peaks(all_peak_af, vcf_d, _bam,  MAX_SNP_IN_CLUSTER = 15, MAX_DEPTH=10000, _out_dir='.', _sample_id='out', _sample_idx='1', is_check_snp_only=True):
    all_rst = []
    max_combined_score = -np.inf
    max_combined_score_i = -1
    for peak_af_i, max_peak_af in enumerate(all_peak_af):
        print('-----------------\nchecking AF %.2f' % (max_peak_af))
#         MAX_SNP_IN_CLUSTER = 10
        _af_dis = 0.05
        _max_af = 0.90
        if max_peak_af == 0:
            _af_dis = 1
        tmp_vcf_d = [i[:] for i in vcf_d]
        if is_check_snp_only:
            tmp_vcf_d = [i for i in tmp_vcf_d if len(i[2]) == 1 and len(i[3]) == 1]
        tmp_vcf_d = [i for i in tmp_vcf_d if abs(i[-1] - max_peak_af) < _af_dis]
        tmp_vcf_d = [i for i in tmp_vcf_d if i[-1] < _max_af]
        
        for _i in range(len(tmp_vcf_d)):
            tmp_vcf_d[_i].append(abs(tmp_vcf_d[_i][-1] - max_peak_af))
        tmp_vcf_d.sort(key=lambda x: x[-1])

        _TMP_MAX_SNP_IN_CLUSTER = min(len(tmp_vcf_d), MAX_SNP_IN_CLUSTER)
        print('candidates %s' % (len(tmp_vcf_d)))
        _HAVE_HIT = False
#         _re_try_t = 1
        _re_try_t = 1 if len(tmp_vcf_d) > 30 else 6
        for tar_vcf_d in _get_sample_times(tmp_vcf_d, _re_try_t, MAX_SNP_IN_CLUSTER):
            if _HAVE_HIT or len(tar_vcf_d) < 5:
                break
            print(len(tar_vcf_d))
            print([(i[0], i[1], i[2], i[3], i[5]) for i in tar_vcf_d])
            print('reading bam from varaints candidates')
            readid_dict, read_cov_dict = get_read_v_matrix(_bam, tar_vcf_d, MAX_DEPTH=MAX_DEPTH)
            print('using %s varaints candidates' % (len(read_cov_dict)))
            _TMP_MAX_SNP_IN_CLUSTER -= 10

            df=pd.DataFrame(readid_dict)
            df.index=df.index.astype(str)
            df=df.sort_index(axis=0).sort_index(axis=1)
            df=df.dropna(axis=1)
            trans=df.T
            trans_list=trans.sum(axis=1).to_list()
            print('get %d read' % (trans.shape[0]))
            trans.head()

            trans_list=trans.sum(axis=1).to_list()
            Z = linkage(trans_list, method='ward')
#             Z = linkage(trans_list, method='centroid', metric='hamming')

            last = Z[-10:, 2]
            last_rev = last[::-1]
            idxs = np.arange(1, len(last) + 1)
            acceleration = np.diff(last, 2)  # 2nd derivative of the distances
    #         tmp_a = np.diff(np.diff(last, 1), 1) / np.diff(last, 1)[1:]
    #         tmp_a[tmp_a < 0] = 0
    #         acceleration =np.diff(last, 2) * tmp_a

            acceleration_rev = acceleration[::-1]
            tar_k = acceleration_rev.argmax() + 2  # if idx 0 is the max of this we want 2 clusters


            fig = plt.figure(figsize = (6, 5))
            ax = fig.add_subplot(1,1,1)
            plt.plot(idxs, last_rev, label='ward\'s distances')
    #         plt.plot(idxs[:-2] + 1, acceleration_rev, label='2nd derivative of the distances')
            plt.plot(idxs[:-2] + 1, acceleration_rev, label='elbow score of the distances')

            ax.scatter(tar_k, acceleration_rev[tar_k-2], c = 'r',marker = 'o', label='selected # of clusters')
            plt.legend()

            plt.xlabel('Number of clusters')
            plt.ylabel('distances')
            sample_id = '%s_%s_AF_%.2f' % (_sample_id, _sample_idx, max_peak_af)
            ax.title.set_text('%s ward\'s distance plot for different clusters' % (sample_id))
    #         plt.savefig('%s/%s.png' % (_out_dir, sample_id))
            plt.savefig('%s/%s_peak_%.2f.png' % (_out_dir, 'distance', max_peak_af))
    #         plt.clf()
    #         print('target k for clustering is %s' % (tar_k))

            _tmp_k = tar_k
            while _tmp_k >= 2:
                print('try using %s clusters' % (_tmp_k))
                n_cluster = _tmp_k
                tag=fcluster(Z, t=n_cluster, criterion='maxclust')

                ori_set = np.array(trans_list)
                all_mean_squared_error = []
                all_centers = []
                for tar_i in range(1, n_cluster+1):
                #     tar_i = 2
                    tar_set = ori_set[np.where(tag == tar_i)]
                    cluster_centers = tar_set.mean(0)
                    all_centers.append(cluster_centers)
                    mse = np.square(tar_set -cluster_centers).mean()
                    all_mean_squared_error.append(mse)
        #             print(tar_i, mse)

                within_c_dis = np.array(all_mean_squared_error).mean()

                all_mean_squared_error_centers = []
                for _i in range(len(all_centers)):
                    for _j in range(1, len(all_centers)):
                        if _j > _i:
                            mse = np.square(all_centers[_i] -all_centers[_j]).mean()
                            all_mean_squared_error_centers.append(mse)
                between_c_dis = np.array(all_mean_squared_error_centers).mean()
                combimed_score = between_c_dis - within_c_dis
                radio_score_ori = between_c_dis * 1. / within_c_dis
                radio_score = radio_score_ori * (len(read_cov_dict)/MAX_SNP_IN_CLUSTER)
                print('generated %d subclusters, intra-clusters dis %.2f, inter-clusters dis %.2f, combined %.2f, radio_ori %.2f, radio %.2f' % (_tmp_k, within_c_dis, between_c_dis, combimed_score, radio_score_ori, radio_score))
#                 if radio_score_ori < 2:
                if radio_score_ori < _tmp_k * 0.75:
                    _tmp_k -= 1
                    continue
                all_rst.append([max_peak_af, tar_vcf_d, readid_dict, read_cov_dict, trans, Z, _tmp_k, combimed_score, radio_score])
                
                _tmp_k -= 1
                _HAVE_HIT = True
                if radio_score > max_combined_score:
                    max_combined_score = radio_score
        #             max_combined_score_i = peak_af_i
                    max_combined_score_i = len(all_rst) - 1
                if radio_score / (_tmp_k + 1) >= 2:
                    print('get high score, exit')
                    return all_rst, max_combined_score_i
    return all_rst, max_combined_score_i


def _get_tar_level_info(_idx, _pct, _tar_level):
    #     _tar_level = 3
        level_k_id = []
        level_k_ptc = []
        for _i in range(len(_idx)):
            _level = int((len(_idx[_i]) - 1) / 2)
            if _level >= _tar_level:
                new_id = _idx[_i][:1+2*_tar_level]
        #         print(_level, new_id, _tar_level)
                if len(level_k_id) > 0 and new_id == level_k_id[-1]:
                    level_k_ptc[-1] += _pct[_i]
                else:
                    level_k_id.append(new_id)
                    level_k_ptc.append(_pct[_i])
            else:
                level_k_id.append('')
                level_k_ptc.append(_pct[_i])
        level_k_ptc = [int(i* 1000) for i in level_k_ptc]
        return level_k_id, level_k_ptc

def generate_clusters_rst_plot(_candidates_list, _out_dir, _sample_id, _top_k):
    _candidates_list.sort(key=lambda x: x[4])
    _idx = [i[-4] for i in _candidates_list]
    _pct = [i[-2] for i in _candidates_list]


    cmap = plt.get_cmap('tab20c')
    colors = cmap(np.arange(100))
#     colors = cmap(np.arange(50))

    fig = plt.figure(figsize = (12, 6))
    ax = fig.add_subplot(1,2,1)


    _i = 1
    _acc_n, _layer_n = 0, 0
    _all_color, _all_label, _all_v = [], [], []
    _lst_label = []
#     while True:
    while _layer_n <= 3:
        labels, v = _get_tar_level_info(_idx, _pct, _i)
        new_label = ['%.0f%%' % (v[i]/10.) if labels[i] != '' else '' for i in range(len(labels))]

        width = 0.3
        wedge_properties = {"width":width, "edgecolor":"w",'linewidth': 2.1}
        patches, texts  = ax.pie(v, labels=new_label, radius=0.3+width*_i, colors = colors[_acc_n:], rotatelabels=False, labeldistance=0.40+0.13*_i, wedgeprops=wedge_properties, startangle=90, counterclock=False)
        [patches[i].set_alpha(0) for i in range(len(patches)) if labels[i] =='']
        _i += 1
        _v_c = sum([1 for i in labels if i != ''])
        if _v_c == 0:
            break
        _all_color = _all_color + [colors[_acc_n+i] for i in range(len(labels)) if labels[i] != '']
        _lst_label = [labels[i] for i in range(len(labels)) if labels[i] != '']
        _all_label = _all_label + [labels[i] for i in range(len(labels)) if labels[i] != '']
        _all_v = _all_v + [v[i] for i in range(len(labels)) if labels[i] != '']

        _acc_n += _v_c
        _layer_n += 1
        
    def check_if_have_suf(x, all_x):
        for row in all_x:
            i = row[4]
            if len(i) > len(x) and x == i[:len(x)]:
                return True
        return False
    
    _all_label_m = [i+'*' if (i in _lst_label and check_if_have_suf(i, _candidates_list)) else i for i in _all_label]

    
    all_batch = [mpatches.Patch(color=_all_color[i], label='Cluster %s' % (_all_label_m[i])) for i in range(_acc_n)]
    plt.legend(handles=all_batch, loc='best', bbox_to_anchor=(-0.3, 0.7, 0.02, 0.2))
    plt.title('%s top clusters compositions' % (_sample_id), y=1.26)


    ax2 = fig.add_subplot(1,2,2)

    # merge all plot infor for subgroup
    p_all_color = [_all_color[i] for i in range(len(_all_label)) if not _all_label[i]+'_1' in _all_label]
    p_all_label = ['%s' % (_all_label_m[i]) for i in range(len(_all_label)) if not _all_label[i]+'_1' in _all_label]
    p_all_v = [_all_v[i]/1000. for i in range(len(_all_label)) if not _all_label[i]+'_1' in _all_label]

    all_info = [list(i) for i in zip(p_all_color, p_all_label, p_all_v)]
    all_info.sort(key=lambda x: -x[2])
    all_info = all_info[:_top_k]
    all_info.append([colors[-1], 'all clusters', sum([i[-1] for i in all_info])])
    ax2.bar([i[1] for i in all_info], [i[2] for i in all_info], tick_label = ['' for i in all_info], color=[i[0] for i in all_info])
    ax2.set_xticks([])
    ax2.set_ylim([0, 1.08])

    _tar_snp_c = []
    for _i in all_info:
        tar_c = [_j[-1][1] for _j in _candidates_list if _j[-4] == _i[1]]
        _tar_snp_c.append('%.0f' % (tar_c[0]) if len(tar_c) > 0 else '-')
    _tar_indel_c = []
    for _i in all_info:
        tar_c = [_j[-1][2] for _j in _candidates_list if _j[-4] == _i[1]]
        _tar_indel_c.append('%.0f' % (tar_c[0]) if len(tar_c) > 0 else '-')
        
    _tar_n_c = []
    for _i in all_info:
        tar_c = [_j[-1][0] for _j in _candidates_list if _j[-4] == _i[1]]
        _tar_n_c.append('%.0f' % (tar_c[0]) if len(tar_c) > 0 else '-')
        
    _tar_af_g9 = []
    for _i in all_info:
        tar_c = ['%.2f' % (_j[-1][3]) for _j in _candidates_list if _j[-4] == _i[1]]
        _tar_af_g9.append(tar_c[0] if len(tar_c) > 0 else '-')
    
    data = [['%.1f' % (i[2]* 100.) for i in all_info],
            _tar_n_c,
            _tar_af_g9,
            _tar_snp_c,
            _tar_indel_c]
    rows = ['Proportions (%)', 'Coverage', 'Median AF', '# of SNP', ' # of INDEL']
    
    columns = [i[1] for i in all_info]
    the_table = plt.table(cellText=data,
                          rowLabels=rows,
                          colLabels=columns,
                          rowLoc='center',
                          cellLoc='center')
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(9)
    the_table.scale(1.0, 1.5)
    # Adjust layout to make room for the table:
    plt.subplots_adjust(bottom=0.3)
    for i, v in enumerate([i[2] for i in all_info]):
        ax2.text(i - .38, v + 0.01, '%.1f%%' % (v*100), fontweight='bold', rotation=45)

    plt.ylabel('proportions')
    plt.title('Final top clusters compositions', y=1.)
    plt.savefig('%s/%s_clustering_rst.png' % (_out_dir, _sample_id), dpi=100, bbox_inches = "tight")
    
# run cluster and return node and their percentage
def run_clustering(_candidate_item, tree, is_tree_empty, _bed, run_clair_path, _n_max_candidates=15, _min_af=0.05, _n_max_coverage=10000, _n_min_supports=50):
    _vcf, _bam, _out_dir, _sample_id, _sample_idx, is_check_clustering, percentage, _v_info = _candidate_item
    print('run', _sample_id, _sample_idx, is_check_clustering, percentage, _v_info)
    os.system(f"mkdir -p {_out_dir}")
    
    if is_check_clustering == 1:
        return [_candidate_item]
    
    _coverage_c, _snp_c, _indel_c, _median_af, _af_1, _af_2, _af_3 = _v_info
    
#     vcf_d, v_cnt = read_vcf(_vcf, is_tree_empty, tree, is_snp_only=True)
    _is_check_snp_only = True
#     _is_check_snp_only = False
    vcf_d, v_cnt = read_vcf(_vcf, is_tree_empty, tree, is_snp_only=_is_check_snp_only)
#     vcf_d = update_af(_bam, _bed, vcf_d)
    _snp_c, _indel_c = v_cnt[0], v_cnt[1]
    all_af = np.array([i[-1] for i in vcf_d])
    af_bins = np.asarray([0, .3, .7, 1.])
    counts, _ = np.histogram(all_af, af_bins)
    _median_af, _af_1, _af_2, _af_3 = np.median(all_af), counts[0], counts[1], counts[2]
    print("read snp, indel, median af, #_af_0_0.3, #_af_0.3_0.7, #_af_0.7_1.0")
    print("%s, %s, %.2f, %d, %d, %d" % ( _snp_c, _indel_c, _median_af, _af_1, _af_2, _af_3))
    
    _v_info = (_coverage_c, _snp_c, _indel_c, _median_af, _af_1, _af_2, _af_3)
    
    if is_check_clustering == 2 or percentage < _min_af or _snp_c < 5:
        return [[_vcf, _bam, _out_dir, _sample_id, _sample_idx, 1, percentage, _v_info]]
    
    rst_peaks, max_peak_af = find_peaks_af(vcf_d, _out_dir, _sample_id, _sample_idx, 0.95, _is_check_snp_only)
    rst_peaks_1, max_peak_af_1 = find_peaks_af(vcf_d, _out_dir, _sample_id, _sample_idx, 0.8, _is_check_snp_only)
    if len(rst_peaks) < 1:
        print('no peaks %s' % (_sample_idx))
        return [[_vcf, _bam, _out_dir, _sample_id, _sample_idx, 1, percentage, _v_info]]
#     print('peaks', rst_peaks)
    all_peak_af = [i[0] for i in rst_peaks]
    all_peak_af_1 = [i[0] for i in rst_peaks_1]
    all_peak_af = list(set(all_peak_af + all_peak_af_1))
    print(all_peak_af)
    all_rst, max_combined_score_i = get_clusers_from_peaks(all_peak_af, vcf_d, _bam, _n_max_candidates, _n_max_coverage, _out_dir, _sample_id, _sample_idx, _is_check_snp_only)
    # clustering fail
    if len(all_rst) < 1:
        print('no clusters %s' % (_sample_idx))
        return [[_vcf, _bam, _out_dir, _sample_id, _sample_idx, 1, percentage, _v_info]]
    max_peak_af, tar_vcf_d, readid_dict, read_cov_dict, trans, Z, tar_k, combimed_score, radio_score = all_rst[max_combined_score_i]
    print('best peak is %.2f, ratio score %.2f' % (max_peak_af, radio_score))
    
    _coverage_c = len(trans)
    _v_info = (_coverage_c, _snp_c, _indel_c, _median_af, _af_1, _af_2, _af_3)
    print('read coverage %s' % (_coverage_c))
    if _coverage_c < _n_min_supports:
        print('no clusters, low coverage %s' % (_sample_idx))
        return [[_vcf, _bam, _out_dir, _sample_id, _sample_idx, 1, percentage, _v_info]]
    trans_s = trans[:]
    trans_list=trans_s.sum(axis=1).to_list()
    Z = linkage(trans_list, 'ward')
    n_cluster = tar_k
    
    tag=fcluster(Z, t=n_cluster, criterion='maxclust')

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(trans_list)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'])
    
    finalDf = pd.concat([principalDf, pd.DataFrame(fcluster(Z, t=n_cluster, criterion='maxclust'))], axis = 1)
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1)
    all_colors=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#000000','#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    colors=all_colors[:n_cluster]
    targets=[i for i in range(1,n_cluster+1)]
    for target, color in zip(targets,colors):
        indicesToKeep = finalDf[0] == target
        scatter = ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1'], finalDf.loc[indicesToKeep, 'principal component 2'], \
                             c = color, s = 1, alpha=0.5, label='Sample #: %s (%.1f%%)' % (len(finalDf.loc[indicesToKeep]), 100.*len(finalDf.loc[indicesToKeep])/len(finalDf)))

    ax.title.set_text('%s, PCA plot' % (_sample_id))
    plt.legend()
    plt.savefig('%s/PCA_AF%.2f.png' % (_out_dir, max_peak_af))
#     plt.clf()
    trans['tag']=tag
    tag_cnt = trans['tag'].value_counts()
#     _min_af = 0.05
#     _tar_f = 0
    _MIN_SUPPORT_R = _n_min_supports
    if str(_sample_idx) != '1':
        for _i in range(1, len(tag_cnt)+1):
            print(tag_cnt.loc[_i] / len(trans) * percentage)
            if tag_cnt.loc[_i] <= _MIN_SUPPORT_R:
                print('subtype limited read support %s, end split' % (tag_cnt.loc[_i]))
                return [[_vcf, _bam, _out_dir, _sample_id, _sample_idx, 1, percentage, _v_info]]
                
#             if tag_cnt.loc[_i] / len(trans) * percentage < _min_af:
#                 print('have clusters proportions < %.3f' % (_min_af))
#                 _tar_f = 2
#                 return [[_vcf, _bam, _out_dir, _sample_id, _sample_idx, 1, percentage, _v_info]]
            
    with open('%s/clusters_info.tsv' % (_out_dir), 'w+') as F:
        F.write('tagged read\t%s\n' % (len(trans)))
        F.write('n of clusters\t%s\n' % (len(tag_cnt)))
        print('genreated %s clusters' % (len(tag_cnt)))
        for _i in range(1, len(tag_cnt)+1):
#             print('cluster %s %s %.0f%%' % (_i, tag_cnt.loc[_i], 100. * tag_cnt.loc[_i] / len(trans)))
            print('cluster %s %s %.2f%%' % (_i, tag_cnt.loc[_i], 100. * tag_cnt.loc[_i] / len(trans) * percentage))
            F.write('cluster_%s coverage\t%s\n' % (_i, tag_cnt.loc[_i]))
    
#     _sample_idx = 1
    trans['tag'].to_csv(f'{_out_dir}/{_sample_idx}_tag.csv',header=False)

    cmd = f"python /autofs/bal31/jhsu/home/projects/HIV/tools/ClusterV/cv/tag_read.py --n_cluster {n_cluster} --out_fn {_out_dir}/{_sample_idx} --split_bam {_bam} {_out_dir}/{_sample_idx}_tag.csv"
    os.system(cmd)
    cmd = 'ls %s/%s.tagged.tag*.bam | parallel "mkdir -p %s/{/}_dir"' % (_out_dir, _sample_idx, _out_dir)
    os.system(cmd)
    
    print('running clair ensenmble')
#     run_clair_path = "/autofs/bal31/jhsu/home/projects/HIV/scripts/run_Clair_ensemble_cv.sh"
    cmd = 'cd %s; ls %s/%s.tagged.tag*.bam |  parallel -j 3 --joblog %s/run_all.log "time bash %s %s/{/} {/} %s/{/}_dir > %s/{/}_dir/run.log"' % \
    (_out_dir, _out_dir, _sample_idx, _out_dir, run_clair_path, _out_dir, _out_dir, _out_dir)
    print(cmd)
    os.system(cmd)
    
    
    new_clusters = []
    for _tag_i in range(1, len(tag_cnt) + 1):
        new_bam = '%s/%s.tagged.tag%s.bam' % (_out_dir, _sample_idx, _tag_i)
        new_vcf = '%s/%s.tagged.tag%s.bam_dir/out.vcf' % (_out_dir, _sample_idx, _tag_i)
        new_sample_idx = '%s_%s' % (_sample_idx, _tag_i)
        new_out_dir = '%s/../%s' % (_out_dir, new_sample_idx)
        new_percentage = 1.0 * tag_cnt.loc[_tag_i] / len(trans) * percentage
        _tar_f = 0
        _tmp_v_info = (tag_cnt.loc[_tag_i], 0, 0, 0, 0, 0, 0)
        new_clusters.append([new_vcf, new_bam, new_out_dir, _sample_id, new_sample_idx, _tar_f, new_percentage, _tmp_v_info])
    print('finish %s, get %s clusters' % (_sample_idx, len(tag_cnt)))
    return new_clusters


def run_CV(args):
    plt.rcParams.update({'figure.max_open_warning': 0})
    random.seed(42)
    
    _vcf = args.vcf_fn
    _bam = args.bam_fn
    _bed = args.bed_fn
    _sample_id = args.sample_id
    _top_k = args.top_k
    _min_af = args.min_af
    _n_max_coverage = args.n_max_coverage
    _n_max_candidates = args.n_max_candidates
    _n_min_supports = args.n_min_supports
    run_clair_path = args.clair_ensemble_s_path
    _out_dir, _bam_n = '/'.join(_vcf.split('/')[:-1]) + '/clustering', _bam.split('/')[-1]
    os.system(f"mkdir -p {_out_dir}")
    
    print(_sample_id, _vcf, _out_dir)
    tree = bed_tree_from(bed_file_path=_bed)
    is_tree_empty = len(tree.keys()) == 0

    # candadates list [vcf, bam, out_dir, is_check_clustering, percentage]
    _candidates_list = []
    _out_dir, _bam_n = '/'.join(_bam.split('/')[:-1]) + '/clustering', _bam.split('/')[-1]
    
#     _coverage_c, _snp_c, _indel_c, _median_af, _af_1, _af_2, _af_3 = 0, 0, 0, 0, 0, 0, 0
#     _v_info = (_coverage_c, _snp_c, _indel_c, _median_af, _af_1, _af_2, _af_3)
    
#     _ori_data = [_vcf, _bam, _out_dir + '/1', _bam_n, 1, 0, 1, _v_info]
#     _candidates_list.append(_ori_data)
#     while len(_candidates_list) < _top_k and any([_i[5] == 0 for _i in _candidates_list]):
#         print('# of generated dataset in candidates: %s' % len(_candidates_list))
# #         print(_candidates_list)
#         # sort by percentage
#         _candidates_list.sort(key=lambda x: -x[-2])
#         _flag_run_cluster = 1
#         _new_candidates_list = []
#         for _k in _candidates_list:
#     #         print(_k)
#             if _k[5] == 0 and _flag_run_cluster == 1:
#                 # run clustering at the current node
#                 # append all sub clusters in to _new_candidates_list
#                 rnt_item = run_clustering(_k, tree, is_tree_empty, _bed, run_clair_path, _n_max_candidates=_n_max_candidates, _min_af=_min_af, _n_max_coverage=_n_max_coverage, _n_min_supports=_n_min_supports)
#                 _new_candidates_list = _new_candidates_list + rnt_item
#     #             print(rnt_item)
#                 _flag_run_cluster = 0
#             else:
#                 _new_candidates_list.append(_k)
#         _candidates_list = _new_candidates_list

#     # update _snp_cnt
#     _new_candidates_list = []
#     for _k in _candidates_list:
#         if _k[5] == 0 or _k[5] == 2:
#             _k[5] = 2
#             rnt_item = run_clustering(_k, tree, is_tree_empty, _bed, run_clair_path)
#             _new_candidates_list = _new_candidates_list + rnt_item
#         else:
#             _new_candidates_list.append(_k)
#     _candidates_list = _new_candidates_list

#     with open('%s/all_clusters_info.tsv' % (_out_dir), 'w+') as F:
#         F.write('%s\n' % ('\t'.join([str(i) for i in ['_vcf', '_bam', '_out_dir', '_sample_id', '_sample_idx', 'is_check_clustering', 'percentage', '_v_info:coverage;snp_c;indel_c;median_af;af1;af2;af3']])))
        
#         for _i in _candidates_list:
#             F.write('%s' % ('\t'.join([str(i) for i in _i[:-1]])))
#             F.write('\t%s\n' % (';'.join([str(i) for i in _i[-1]])))
            
    _candidates_list = []
    with open('%s/all_clusters_info.tsv' % (_out_dir), 'r') as F:
        for line in F:
            if len(line) > 0 and line[0] != '_':
                row = line.split()
                n_row = row[:6] + [float(row[6])] + row[7:-1] + [[float(i) for i in row[-1].split(';')]]
#                 print(n_row)
                _candidates_list.append(n_row)        
    print(_candidates_list)
    
    if len(_candidates_list) > 1:
        generate_clusters_rst_plot(_candidates_list, _out_dir, _sample_id, _top_k)

    return _candidates_list


def main():
    parser = ArgumentParser(description="Cluster alignment based on variant")

    parser.add_argument('--vcf_fn', type=str, default="input.vcf",
                    help="Truth VCF file input, required")

    parser.add_argument('--bam_fn', type=str, default="input.bam",
                        help="input bam file, required")

    parser.add_argument('--bed_fn', type=str, default="input.bed",
                        help="input bed file, optional")
    
    parser.add_argument('--sample_id', type=str, default="sample_id",
                        help="sample_id, optional")
    
    parser.add_argument('--top_k', type=int, default=50,
                        help="top_k, optional")
    
    parser.add_argument('--min_af', type=float, default=0.025,
                        help="min_af, optional")
    
    parser.add_argument('--n_max_coverage', type=int, default=10000,
                        help="max read for clustering, optional")
    
    parser.add_argument('--n_max_candidates', type=int, default=15,
                        help="max selected candidates for clustering, optional")
    
    parser.add_argument('--n_min_supports', type=int, default=30,
                        help="minimum read support for creating a subtype, optional")

    parser.add_argument('--clair_ensemble_s_path', type=str, default="/autofs/bal31/jhsu/home/projects/HIV/scripts/run_Clair_ensemble_cv.sh",
                        help="clair ensemble scipt path")
    
    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    run_CV(args)

if __name__ == "__main__":
    main()
