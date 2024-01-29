# -*- coding: UTF-8 -*
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier


call_name_dic = {"bcftools": 1, "Delly": 2, "freebayes": 3, "gatk": 4, "Manta": 5, "Pindel": 6, "SvABA": 7}
snp_caller_arr = ["bcftools", "freebayes", "gatk"]
sv_caller_arr = ["Delly", "Manta", "Pindel", "SvABA"]
snp_caller_sets = [["bcftools", "freebayes"], ["freebayes", "gatk"], ["bcftools", "gatk"], ["bcftools", "freebayes", "gatk"]]
sv_caller_sets = [["Delly", "Manta"], ["Delly", "Pindel"], ["Delly", "SvABA"], ["Manta", "Pindel"], ["Manta", "SvABA"], ["Pindel", "SvABA"], ["Delly", "Manta", "Pindel"], ["Delly", "Manta", "SvABA"], ["Delly", "Pindel", "SvABA"], ["Manta", "Pindel", "SvABA"], ["Delly", "Manta", "Pindel", "SvABA"]]

def read_file(file_dir):
    readFileObject = open(file_dir)
    try:
        readFile = readFileObject.read( )
    finally:
        readFileObject.close( )
    readFile_rows = readFile.split('\n')
    readFile_rows.pop()
    return readFile_rows

def get_feature(feature_file_rows):
  all_row_feature_file_rows_arr = []
  for one_feature_file_row in feature_file_rows:
    ini_feature_arr = one_feature_file_row.split(",split,")
    str_feature_arr = ini_feature_arr[0].split(",")[2:]
    feature_arr = []
    for one_str_feature_arr_ele in str_feature_arr:
        feature_arr.append(float(one_str_feature_arr_ele))
    all_row_feature_file_rows_arr.append(feature_arr)
  return all_row_feature_file_rows_arr

def get_target(target_file_rows):
  all_row_snp_caller_target_arr = []
  all_row_sv_caller_target_arr = []
  for one_target_file_row in target_file_rows:
    snp_caller_target_arr = []
    sv_caller_target_arr = []
    this_sample_allCaller_goodToBad_arr = one_target_file_row.split(",")[1:]
    this_sample_snpCaller_goodToBad_arr = []
    this_sample_svCaller_goodToBad_arr = []
    for one_this_sample_allCaller_goodToBad_arr_ele in this_sample_allCaller_goodToBad_arr:
        if one_this_sample_allCaller_goodToBad_arr_ele in snp_caller_arr:
            this_sample_snpCaller_goodToBad_arr.append(one_this_sample_allCaller_goodToBad_arr_ele)
        if one_this_sample_allCaller_goodToBad_arr_ele in sv_caller_arr:
            this_sample_svCaller_goodToBad_arr.append(one_this_sample_allCaller_goodToBad_arr_ele)
    #print this_sample_snpCaller_goodToBad_arr, this_sample_svCaller_goodToBad_arr
    for one_snp_caller_set in snp_caller_sets:
        this_set_index_arr = []
        for one_snp_caller in one_snp_caller_set:
            this_set_index_arr.append(this_sample_snpCaller_goodToBad_arr.index(one_snp_caller))
        this_dic = dict(zip(this_set_index_arr, one_snp_caller_set))
        min_index = min(this_set_index_arr)
        this_snp_target = this_dic[min_index]
        snp_caller_target_arr.append(this_snp_target)
    for one_sv_caller_set in sv_caller_sets:
        this_set_index_arr = []
        for one_sv_caller in one_sv_caller_set:
            this_set_index_arr.append(this_sample_svCaller_goodToBad_arr.index(one_sv_caller))
        this_dic = dict(zip(this_set_index_arr, one_sv_caller_set))
        min_index = min(this_set_index_arr)
        this_sv_target = this_dic[min_index]
        sv_caller_target_arr.append(this_sv_target)
    #return snp_caller_target_arr, sv_caller_target_arr
    all_row_snp_caller_target_arr.append(snp_caller_target_arr)
    all_row_sv_caller_target_arr.append(sv_caller_target_arr)
  return all_row_snp_caller_target_arr, all_row_sv_caller_target_arr

snp_feature_file_rows = read_file("mData_snp_allWindowMetaFeature.csv")
sv_feature_file_rows = read_file("mData_sv_allWindowMetaFeature.csv")
FPR_target_file_rows = read_file("mData_allCaller_FPR_goodToBad.csv")
FNR_target_file_rows = read_file("mData_allCaller_FNR_goodToBad.csv")

snp_feature_arr = get_feature(snp_feature_file_rows)
sv_feature_arr = get_feature(sv_feature_file_rows)
FPR_snp_target_arr, FPR_sv_target_arr =  get_target(FPR_target_file_rows)
FNR_snp_target_arr, FNR_sv_target_arr =  get_target(FNR_target_file_rows)
'''
one_arr = []
for oneone in FPR_snp_target_arr:
    one_arr.append(oneone[-1])
print one_arr
print len(one_arr)
print snp_feature_arr[0]


clf = RandomForestClassifier(n_estimators=70, max_depth=7, min_samples_split=2, min_samples_leaf=1, max_features=5, oob_score=False, random_state=10)
clf.fit(snp_feature_arr, one_arr)
joblib.dump(clf, "tttttttyyyyyyy.m")
X_test = [[64.7028238352,27.6467903421,38.0600190909,26.4398406375,1.24887102364e-06,0.000311502219453,0.000335203475281,0.0610171901823]]
y_pred = clf.predict(X_test)
print y_pred
'''
clf = RandomForestClassifier(n_estimators=70, max_depth=7, min_samples_split=2, min_samples_leaf=1, max_features=5, oob_score=False, random_state=10)
def getModels(feature_arr, target_arr, this_name_prefix):
    if "sv" in this_name_prefix:
        caller_sets = sv_caller_sets
    if "snp" in this_name_prefix:
        caller_sets = snp_caller_sets
    for i in range(len(caller_sets)):
        this_target = []
        for one_target in target_arr:
            this_target.append(one_target[i])
        clf.fit(feature_arr, this_target)
        this_name_arr = caller_sets[i]
        this_name = ""
        for one_this_name_arr_ele in this_name_arr:
            this_name = this_name + str(call_name_dic[one_this_name_arr_ele])
        joblib.dump(clf, this_name_prefix + this_name + ".m")
getModels(snp_feature_arr, FPR_snp_target_arr, "FPRsnp_")
getModels(snp_feature_arr, FNR_snp_target_arr, "FNRsnp_")
getModels(sv_feature_arr, FPR_sv_target_arr, "FPRsv_")
getModels(sv_feature_arr, FNR_sv_target_arr, "FNRsv_")


