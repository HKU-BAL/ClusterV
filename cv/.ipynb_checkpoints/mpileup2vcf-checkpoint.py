import argparse
import sys
import os

# a simple script to get one read's mpileup and transfer it from fasta to vcf
# e.g samtools mpileup -f ref.fasta bam | python mpileup2vcf.py | bcftools view -O z -o test.vcf.gz

def _get_header(_sample_n):
    _s = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Confidence in this variant being real is below calling threshold.">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=LENGUESS,Number=.,Type=Integer,Description="Best guess of the indel length">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1)">
##contig=<ID=NC_001802.1,length=9181>
##bcftools_filterVersion=1.9+htslib-1.9
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	%s
""" % (_sample_n)
    return _s

def _get_alt_indel(_s):
    base_idx = 1
    advance = 0
    while True:
        num = _s[base_idx]
        if num.isdigit():
            advance = advance * 10 + int(num)
            base_idx += 1
        else:                                                                                                                                                                                                                               break
    _n, _alt = advance, _s[base_idx:base_idx+advance]
    return _n, _alt

def main(args):
    
    all_rec = []
    for line in args.input_f:
        line = line.strip()
        if len(line) <= 1:
            continue
        row = line.split('\t')
        _contig, _pos, _ref, _o_alt = row[0], int(row[1]), row[2], row[4].strip()
        _tmp_alt = [k for k in _o_alt.split(',') if len(k) > 0]    
        _alt = '' if len(_tmp_alt)< 1 else _tmp_alt[0]
            
        if len(_alt) <= 0 or _alt in ['.', '*']:
            continue
        if _alt[0] == '^' or _alt[-1] == '$':
            continue
#         print("%s %s" % (_o_alt.split(','), _tmp_alt), file=sys.stderr)
#         print(line, file=sys.stderr)
        if (_alt[0] == '.' and _alt[1] in ['-', '+']) or (_alt[0] in ['-', '+']):
            if _alt[0] == '.':
                _alt = _alt[1:]
            # ins
            _n, _n_alt = _get_alt_indel(_alt)
#             print(_n, _n_alt)
            if _alt[0] == '+':
                _alt = _ref + _n_alt
            else:
                _alt = _ref
                _ref = _ref + _n_alt
        elif len(_alt) > 1:
            print('warning at site %s' % row, file=sys.stderr)
            _alt = _alt[0]
        all_rec.append([_contig, _pos, _ref, _alt])
#         print(_contig, _pos, _ref, _alt)
    print(_get_header(args.sample_n), end='')
    for row in all_rec:
        _info_row = [row[0], str(row[1]), '.', row[2].upper(), row[3].upper(), '100', 'PASS', '.', 'GT:GQ:DP:AF', '1/1:100:100:1']
        print("%s" % ('\t'.join(_info_row)))
        


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--out_f', type=str, default='out.vcf')
    parser.add_argument('--input_f', type=argparse.FileType('r'), default=(None if sys.stdin.isatty() else sys.stdin))
    parser.add_argument('--sample_n', type=str, default='1')

    args = parser.parse_args()
    main(args)
