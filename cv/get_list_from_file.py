def get_list_from_file(file):
    with open(file) as fi:
        list=[line.strip() for line in fi]
    return list
def get_list_from_csv(file):
    with open(file) as fi:
        list=[line.strip().split(',') for line in fi]
    return list
def get_dict_from_csv(file):
    d={}
    with open(file) as fi:
        for line in fi:
            key,value=line.strip().split(',') 
            d[key]=value
    return d
def get_set_from_file(file):
    with open(file) as fi:
        set={line.strip() for line in fi}
    return set

