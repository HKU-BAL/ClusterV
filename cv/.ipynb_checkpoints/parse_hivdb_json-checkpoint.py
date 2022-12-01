import json

# parse json, reture a list of ['gene', 'drugClass', 'drugName', 'drugScore', 'resistanceLevel', 'mutation', 'mutationScore', 'mutationType', 'comments']
def parse_hivdb_json(_data):
    all_l = []
    all_l.append(['gene', 'drugClass', 'drugName', 'drugScore', 'resistanceLevel', 'mutation', 'mutationScore', 'mutationType', 'comments'])
    _tar_json = None
    try:
        _tar_json = _data[0]
    except:
        _tar_json = _data
        
    if 'drugResistance' in _tar_json:
        for i in _tar_json['drugResistance']:
        #     print(i['gene']['name'])
            gene_name = i['gene']['name']
            for j in i['drugScores']:
                _drug_n = j['drug']['name']
                _drug_cls = j['drugClass']['name']
                _drug_score = j['score']
                _drug_resistance_level = j['text']

                # mutation, mutationType, comment
                mutations_l = []
                for ps in j['partialScores']:
                    p_score = ps['score']
                    for mut in ps['mutations']:
                        mutations_l.append([mut['text'], mut['primaryType'], mut['comments'][0]['text']])
                        # gene name, drug cls, drug name, drug score, mutation score, mutation text, mutation type, mutations comment
    #                     _rst_l = [gene_name, _drug_cls, _drug_n, _drug_score, mut['text'], p_score, mut['primaryType']]
                        _cmt = mut['comments'][0]['text'].strip().replace("\t", ";")
                        _rst_l = [gene_name, _drug_cls, _drug_n, _drug_score, _drug_resistance_level, mut['text'], p_score, mut['primaryType'], _cmt]
                        all_l.append(_rst_l)
    # get other variants
    try: 
        all_v_m = [i[5] for i in all_l[1:]]
        for i in _tar_json['alignedGeneSequences']:
            gene_name = i['gene']['name']
            for j in i['mutations']:
                _mutation = j['text']
                _type = j['primaryType']
                if _mutation not in all_v_m:
                    _rst_l = [gene_name, 'NA', 'NA', 'NA', 'NA', _mutation, 'NA', _type, 'NA']
                    all_l.append(_rst_l)
    except:
        pass
    return all_l

# with open(_out_f, 'r') as f:
#     _data = json.load(f)
#     all_l = parse_hivdb_json(_data)
#     for i in all_l:
#         print(i)