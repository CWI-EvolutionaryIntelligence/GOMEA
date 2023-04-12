import json
            
def writeResultsToJson(input,prob_tag,lang_tag,lm_tag):
    res_dir = "results/"
    if prob_tag == "trap5":
        res_dir = "deceptive_trap/" + res_dir
    elif prob_tag == "maxcut":
        res_dir = "maxcut/" + res_dir
    full_tag = prob_tag+"_"+lang_tag+"_"+lm_tag
    filename = res_dir+full_tag+".json"
    with open(filename,'w') as fp:
        json.dump(input,fp)

def readResultsFromJson(prob_tag,lang_tag,lm_tag):
    res_dir = "results/"
    if prob_tag == "trap5":
        res_dir = "deceptive_trap/" + res_dir
    elif prob_tag == "maxcut":
        res_dir = "maxcut/" + res_dir
    full_tag = prob_tag+"_"+lang_tag+"_"+lm_tag
    filename = res_dir+full_tag+".json"
    try:
        with open(filename,'r') as fp:
            res = json.load(fp)
            for k,sub_dict in res.items():
                try:
                    res[k] = {int(dim):[float(r) for r in res] for dim,res in sub_dict.items()}
                except Exception as e:
                    res[k] = {int(dim):res for dim,res in sub_dict.items()}
            return res
    except Exception as e:
        print(e)
        return {}

def resultsExist(prob_tag,lang_tag,lm_tag):
    res_dir = "results/"
    if prob_tag == "trap5":
        res_dir = "deceptive_trap/" + res_dir
    elif prob_tag == "maxcut":
        res_dir = "maxcut/" + res_dir
    full_tag = prob_tag+"_"+lang_tag+"_"+lm_tag
    filename = res_dir+full_tag+".json"
    try:
        with open(filename,'r') as fp:
            return True
    except FileNotFoundError:
        return False