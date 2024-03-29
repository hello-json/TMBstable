# -*- coding: UTF-8 -*-
import os
import pysam
import random
import csv
import argparse
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import math

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam_file', help='input the bam file, eg: sample.bam', required=True, type=str)
parser.add_argument('-d', '--ini_info_dir', help='input the ini_info_dir, eg: ini_info/', required=True, type=str)
parser.add_argument('-g', '--refGenome', help='input the refGenome file, eg: refGenome/hg19.fa', required=True, type=str)
parser.add_argument('-m', '--snp_metaModel', help='the snp meta model file dir, eg: Nsnp_134.m', required=True, type=str)
parser.add_argument('-n', '--sv_metaModel', help='the sv meta model file dir, eg: Nsv_2567.m', required=True, type=str)
parser.add_argument('-w', '--snp_window', help='the snp window size, eg: 1000000', required=True, type=str)
parser.add_argument('-r', '--iniRepeatFileDir', help='the repeat file dir, eg: hg19RepeatFile/rmsk.txt', required=True, type=str)
parser.add_argument('-o', '--outFileDir', help='the out file dir, eg: myresult/TMBstable.vcf', required=True, type=str)
args = parser.parse_args()
bamFile = args.bam_file
ini_info_dir = args.ini_info_dir
refGenome = args.refGenome
snp_metaModel = args.snp_metaModel
sv_metaModel = args.sv_metaModel
snp_window = args.snp_window
iniRepeatFileDir = args.iniRepeatFileDir
outFileDir = args.outFileDir


# begin to load model
snp_model = joblib.load(snp_metaModel)
sv_model = joblib.load(sv_metaModel)
# end to load model

chromosomesCodingDic = {"chr1": 1000000000, "chr2": 2000000000, "chr3": 3000000000, "chr4": 4000000000, "chr5": 5000000000, "chr6": 6000000000, "chr7": 7000000000, "chr8": 8000000000, "chr9": 9000000000, "chr10": 10000000000, "chr11": 11000000000, "chr12": 12000000000, "chr13": 13000000000, "chr14": 14000000000, "chr15": 15000000000, "chr16": 16000000000, "chr17": 17000000000, "chr18": 18000000000, "chr19": 19000000000, "chr20": 20000000000, "chr21": 21000000000, "chr22": 22000000000, "chrX": 23000000000, "chrY": 24000000000, "1": 1000000000, "2": 2000000000, "3": 3000000000, "4": 4000000000, "5": 5000000000, "6": 6000000000, "7": 7000000000, "8": 8000000000, "9": 9000000000, "10": 10000000000, "11": 11000000000, "12": 12000000000, "13": 13000000000, "14": 14000000000, "15": 15000000000, "16": 16000000000, "17": 17000000000, "18": 18000000000, "19": 19000000000, "20": 20000000000, "21": 21000000000, "22": 22000000000, "X": 23000000000, "Y": 24000000000}


snp_window_arr = []
get_bed_cmd = "bedtools bamtobed -i " + bamFile + " > tmp.bed"
os.system(get_bed_cmd)
iniBedObject = open("tmp.bed")
try:
     iniBedFile = iniBedObject.read( )
finally:
     iniBedObject.close( )
iniBedFile_rows = iniBedFile.split('\n')
iniBedFile_rows.pop()
window_size_add = 0
for one_iniBedFile_row in iniBedFile_rows:
    one_iniBedFile_row_arr = one_iniBedFile_row.split('\t')
    if one_iniBedFile_row_arr[0] in chromosomesCodingDic:
        one_iniBedFile_row_begin = chromosomesCodingDic[one_iniBedFile_row_arr[0]] + int(one_iniBedFile_row_arr[1])
        one_iniBedFile_row_end = chromosomesCodingDic[one_iniBedFile_row_arr[0]] + int(one_iniBedFile_row_arr[2])
    else:
        continue
    if window_size_add == 0:
        this_window_begin = one_iniBedFile_row_begin
        this_window_end = one_iniBedFile_row_end
        window_size_add = window_size_add + (one_iniBedFile_row_end - this_window_begin)
        store_this = [this_window_begin, this_window_end]
        if window_size_add > snp_window:
            window_size_add = 0
            snp_window_arr.append([this_window_begin, this_window_end])
    else:
        if int(str(this_window_end)[0:-9] + "000000000") != chromosomesCodingDic[one_iniBedFile_row_arr[0]]:
            snp_window_arr.append([this_window_begin, this_window_end])
            window_size_add = 0
            continue
        window_size_add = window_size_add + (one_iniBedFile_row_end - max(this_window_end, one_iniBedFile_row_begin))
        this_window_end = max(this_window_end, one_iniBedFile_row_end)
        if window_size_add > snp_window:
            snp_window_arr.append([this_window_begin, this_window_end])
            window_size_add = 0
if len(snp_window_arr) == 0:
    snp_window_arr.append(store_this)
window_readLen_arr = []
store_for_mappingAndBaseQual_step_arr = []
store_for_mappingAndBaseQual_step_oneDimensional_arr = []
for one_snp_window in snp_window_arr:
    this_window_read_count = 0
    this_window_readLen_sum = 0
    this_window_store_for_mappingAndBaseQual_arr = []
    for one_iniBedFile_row in iniBedFile_rows:
        one_iniBedFile_row_arr = one_iniBedFile_row.split('\t')
        if one_iniBedFile_row_arr[0] in chromosomesCodingDic:
            one_iniBedFile_row_begin = chromosomesCodingDic[one_iniBedFile_row_arr[0]] + int(one_iniBedFile_row_arr[1])
            one_iniBedFile_row_end = chromosomesCodingDic[one_iniBedFile_row_arr[0]] + int(one_iniBedFile_row_arr[2])
        else:
            continue
        if one_iniBedFile_row_begin > one_snp_window[1]:
            break
        if one_iniBedFile_row_begin > one_snp_window[0] and one_iniBedFile_row_end < one_snp_window[1]:
            this_window_readLen_sum = this_window_readLen_sum + (one_iniBedFile_row_end - one_iniBedFile_row_begin)
            this_window_read_count = this_window_read_count + 1
            if this_window_read_count == 1:
                this_window_store_for_mappingAndBaseQual_arr.append(one_iniBedFile_row_arr[3].split('/')[0])
                store_for_mappingAndBaseQual_step_oneDimensional_arr.append(one_iniBedFile_row_arr[3].split('/')[0])
            else:
                if random.random() > 0.99:        
                    this_window_store_for_mappingAndBaseQual_arr.append(one_iniBedFile_row_arr[3].split('/')[0])
                    store_for_mappingAndBaseQual_step_oneDimensional_arr.append(one_iniBedFile_row_arr[3].split('/')[0])
    try:
        window_readLen_arr.append(float(this_window_readLen_sum) / float(this_window_read_count))
    except:
        window_readLen_arr.append(float(50))
    store_for_mappingAndBaseQual_step_arr.append(this_window_store_for_mappingAndBaseQual_arr)
window_mappingQuality_arr = []
window_sequencingQuality_arr = []
store_for_mappingAndBaseQual_step_oneDimensional_arr_info_dic = {}
sam_rd = pysam.AlignmentFile(bamFile, 'rb')
for rec in sam_rd:
    if rec.query_name in store_for_mappingAndBaseQual_step_oneDimensional_arr:    
        store_for_mappingAndBaseQual_step_oneDimensional_arr_info_dic[rec.query_name] = [rec.mapping_quality, float(sum(rec.query_qualities)) / float(len(rec.query_qualities))]
sam_rd.close()
for one_store_for_mappingAndBaseQual_step_arr_element in store_for_mappingAndBaseQual_step_arr:
  try:
    window_Quality_count = 0
    window_mappingQuality_sum = 0
    window_sequencingQuality_sum = 0
    for one_store_for_mappingAndBaseQual_step_arr_element_element in one_store_for_mappingAndBaseQual_step_arr_element:
        window_Quality_count = window_Quality_count + 1
        window_mappingQuality_sum = window_mappingQuality_sum + store_for_mappingAndBaseQual_step_oneDimensional_arr_info_dic[one_store_for_mappingAndBaseQual_step_arr_element_element][0]
        window_sequencingQuality_sum = window_sequencingQuality_sum + store_for_mappingAndBaseQual_step_oneDimensional_arr_info_dic[one_store_for_mappingAndBaseQual_step_arr_element_element][1]
    window_mappingQuality_arr.append(float(window_mappingQuality_sum) / float(window_Quality_count))
    window_sequencingQuality_arr.append(float(window_sequencingQuality_sum) / float(window_Quality_count))
  except:
    window_mappingQuality_arr.append(float(38))
    window_sequencingQuality_arr.append(float(29))
window_depth_arr = []
window_insRatio_arr = []
window_delRatio_arr = []
window_subRatio_arr = []
pileup_tsv_dic = {"A": 4, "C": 5, "G": 6, "T": 7, "a": 4, "c": 5, "g": 6, "t": 7, "N": 4, "n": 4}
for one_snp_window in snp_window_arr:
  try:
    this_window_begin = one_snp_window[0]
    this_window_end = one_snp_window[1]
    this_chr_num_str = str(this_window_end)[0:-9]
    if this_chr_num_str == "23":
        this_chr = "chrX"
        #this_chr = "X"      #other genome
    elif this_chr_num_str == "24":
        this_chr = "chrY"
        #this_chr = "Y"     #other genome
    else:
        this_chr = "chr" + this_chr_num_str
        #this_chr = this_chr_num_str        #other genome
    this_begin = str(this_window_begin)[-9:]
    this_end = str(this_window_end)[-9:]
    os.system("quickgt " + bamFile + " " + refGenome + " " + this_chr + ":" + this_begin + "-" + this_end + " > tmp.tsv")
    depth_sum = 0
    count = 0
    ins_count = 0
    del_count = 0
    sub_count = 0
    with open('tmp.tsv') as f:
        tsvreader = csv.reader(f, delimiter = '\t')
        for line in tsvreader:
            if line[0] == 'CHR':
                continue
            count = count + 1
            depth_sum = depth_sum + int(line[3])
            if line[9] != "0;0":
                ins_count = ins_count + 1
            elif line[10] != "0;0":
                del_count = del_count + 1
            elif int(line[3]) > (int(line[pileup_tsv_dic[line[2]]].split(";")[0]) + int(line[pileup_tsv_dic[line[2]]].split(";")[1])):
                sub_count = sub_count + 1
    depth = float(depth_sum) / float(count)
    insRatio = float(ins_count) / float(count)
    delRatio = float(del_count) / float(count)
    subRatio = float(sub_count) / float(count)
    window_depth_arr.append(depth)
    window_insRatio_arr.append(insRatio)
    window_delRatio_arr.append(delRatio)
    window_subRatio_arr.append(subRatio)
  except:
    window_depth_arr.append(float(30))
    window_insRatio_arr.append(float(0.0017))
    window_delRatio_arr.append(float(0.0012))
    window_subRatio_arr.append(float(0.0336))
iniRepeatObject = open(iniRepeatFileDir)
try:
     iniRepeatFile = iniRepeatObject.read( )
finally:
     iniRepeatObject.close( )
iniRepeatFile_rows = iniRepeatFile.split('\n')
iniRepeatFile_rows.pop()
repeat_arr = []
for one_iniRepeatFile_row in iniRepeatFile_rows:
    one_iniRepeatFile_row_arr = one_iniRepeatFile_row.split('\t')
    try:                          
        repeat_arr.append([chromosomesCodingDic[one_iniRepeatFile_row_arr[5]] + int(one_iniRepeatFile_row_arr[6]), chromosomesCodingDic[one_iniRepeatFile_row_arr[5]] + int(one_iniRepeatFile_row_arr[7])])
    except:
        continue
window_repeatRatio_arr = []
for one_window in snp_window_arr:
    windowInRepeatLengthRatio = float(0)
    for one_repeat in repeat_arr:
        if one_window[0] >= one_repeat[0] and one_window[0] <= one_repeat[1]:
            windowInRepeatLengthRatio = float(min(one_window[1], one_repeat[1]) - one_repeat[0]) / float(one_window[1] - one_window[0])
            break
    window_repeatRatio_arr.append(windowInRepeatLengthRatio)


sv_window_span = int(snp_window / 2)
all_chr_region_dic = {}
for one_snp_window in snp_window_arr:
    if str(one_snp_window[0])[0:-9] not in all_chr_region_dic:
        all_chr_region_dic[str(one_snp_window[0])[0:-9]] = [one_snp_window]
    else:
        all_chr_region_dic[str(one_snp_window[0])[0:-9]].append(one_snp_window)
for one_chr_region_dic_key in all_chr_region_dic:
    if len(all_chr_region_dic[one_chr_region_dic_key]) > 1:
        all_chr_region_dic[one_chr_region_dic_key] = [all_chr_region_dic[one_chr_region_dic_key][0][0], all_chr_region_dic[one_chr_region_dic_key][-1][-1]]

chr_flag_arr = []
sam_rd = pysam.AlignmentFile(bamFile, 'rb')
for rec in sam_rd:
  if rec.reference_name != None:
    if "chr" in rec.reference_name:
        rec_reference_name = rec.reference_name[3:]
    else:
        rec_reference_name = rec.reference_name
    if rec_reference_name in all_chr_region_dic and rec.cigarstring != None and 'S' in rec.cigarstring and rec_reference_name not in chr_flag_arr:
        rec_pos = int(rec_reference_name) * 1000000000 + rec.pos
        if (rec_pos + sv_window_span) > all_chr_region_dic[rec_reference_name][0][0] and (rec_pos + sv_window_span) < all_chr_region_dic[rec_reference_name][0][1]:
            this_last = all_chr_region_dic[rec_reference_name][1:]
            all_chr_region_dic[rec_reference_name] = [[all_chr_region_dic[rec_reference_name][0][0], (rec_pos + sv_window_span)], [(rec_pos + sv_window_span+1), all_chr_region_dic[rec_reference_name][0][1]]] + this_last
        chr_flag_arr.append(rec_reference_name)
sam_rd.close()

sv_window_arr = []
for i in range(1, 25):
    if str(i) in all_chr_region_dic:
        sv_window_arr = sv_window_arr + all_chr_region_dic[str(i)]
window_length_arr = []
for one_sv_window in sv_window_arr:
    window_length_arr.append(int(one_sv_window[1]) - int(one_sv_window[0]))
window_insertSize_arr = []
for one_sv_window in sv_window_arr:
  try:
    this_sv_window_chr = "chr" + str(one_sv_window[0])[0:-9]
    #this_sv_window_chr = str(one_sv_window[0])[0:-9]    #other genome
    if this_sv_window_chr == "chr23":
        this_sv_window_chr = "chrX"
    if this_sv_window_chr == "23":
        this_sv_window_chr = "X"
    if this_sv_window_chr == "chr24":
        this_sv_window_chr = "chrY"
    if this_sv_window_chr == "24":
        this_sv_window_chr = "Y"
    this_sv_window_begin = str(one_sv_window[0])[-9:]
    this_sv_window_end = str(one_sv_window[1])[-9:]
    os.system("samtools view " + bamFile + " " + this_sv_window_chr + ":" + this_sv_window_begin + "-" + this_sv_window_end + " | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > tmp_insertSize.txt")
    if os.path.getsize("tmp_insertSize.txt") != 0:
        insertSizeFileDir = "tmp_insertSize.txt"
        insertSizeObject = open(insertSizeFileDir)
        try:
            insertSizeFile = insertSizeObject.read( )
        finally:
            insertSizeObject.close( )
        insertSizeFile_rows = insertSizeFile.split('\n')
        insertSizeFile_rows.pop()
        insertSize_sum = 0
        insertSize_account = 0
        for one_insertSizeFile_row in insertSizeFile_rows:
            one_insertSizeFile_row_arr = one_insertSizeFile_row.split(" ")
            insertSize_sum = insertSize_sum + int(one_insertSizeFile_row_arr[0])
            insertSize_account = insertSize_account + 1
        this_window_insertSize_value = float(insertSize_sum) / float(insertSize_account)
    else:
        this_window_insertSize_value = float(0.1)
    window_insertSize_arr.append(this_window_insertSize_value)
  except:
    window_insertSize_arr.append(float(0.1))
window_SVrepeatRatio_arr = []
for one_window in sv_window_arr:
    windowInRepeatLengthRatio = float(0)
    for one_repeat in repeat_arr:         
        if one_window[0] >= one_repeat[0] and one_window[0] <= one_repeat[1]:
            windowInRepeatLengthRatio = float(min(one_window[1], one_repeat[1]) - one_repeat[0]) / float(one_window[1] - one_window[0])
            if (one_window[1] - one_window[0]) < 5500:          
                windowInRepeatLengthRatio = min(float(1), windowInRepeatLengthRatio * 1.5)
            break
    window_SVrepeatRatio_arr.append(windowInRepeatLengthRatio)
window_SVreadLen_arr = []
for one_snp_window in sv_window_arr:
    this_window_read_count = 0
    this_window_readLen_sum = 0
    this_window_store_for_mappingAndBaseQual_arr = []
    for one_iniBedFile_row in iniBedFile_rows:
        one_iniBedFile_row_arr = one_iniBedFile_row.split('\t')
        if one_iniBedFile_row_arr[0] in chromosomesCodingDic:
            one_iniBedFile_row_begin = chromosomesCodingDic[one_iniBedFile_row_arr[0]] + int(one_iniBedFile_row_arr[1])
            one_iniBedFile_row_end = chromosomesCodingDic[one_iniBedFile_row_arr[0]] + int(one_iniBedFile_row_arr[2])
        else:
            continue
        if one_iniBedFile_row_begin > one_snp_window[1]:
            break
        if one_iniBedFile_row_begin > one_snp_window[0] and one_iniBedFile_row_end < one_snp_window[1]:
            this_window_readLen_sum = this_window_readLen_sum + (one_iniBedFile_row_end - one_iniBedFile_row_begin)
            this_window_read_count = this_window_read_count + 1
    try:
        window_SVreadLen_arr.append(float(this_window_readLen_sum) / float(this_window_read_count))
    except:
        window_SVreadLen_arr.append(-1)

missingValueFilling_sum = 0
missingValueFilling_count = 0
for one_window_SVreadLen_arr_value in window_SVreadLen_arr:
    if one_window_SVreadLen_arr_value != -1:
        missingValueFilling_sum = missingValueFilling_sum + one_window_SVreadLen_arr_value
        missingValueFilling_count = missingValueFilling_count + 1
try:
    missingValueFilling = float(missingValueFilling_sum) / float(missingValueFilling_count)
except:
    missingValueFilling = float(10)
tmp_window_SVreadLen_arr = window_SVreadLen_arr
window_SVreadLen_arr = []
for one_tmp_window_SVreadLen_arr_element in tmp_window_SVreadLen_arr:
    if one_tmp_window_SVreadLen_arr_element == -1:
        one_tmp_window_SVreadLen_arr_element = missingValueFilling
    window_SVreadLen_arr.append(one_tmp_window_SVreadLen_arr_element)
window_SVdepth_arr = []
window_SVinsRatio_arr = []
window_SVdelRatio_arr = []
window_SVsubRatio_arr = []
for one_snp_window in sv_window_arr:
    this_window_begin = one_snp_window[0]
    this_window_end = one_snp_window[1]
    this_chr_num_str = str(this_window_end)[0:-9]
    if this_chr_num_str == "23":
        this_chr = "chrX"
        #this_chr = "X"    #other genome
    elif this_chr_num_str == "24":
        this_chr = "chrY"
        #this_chr = "Y"     #other genome
    else:
        this_chr = "chr" + this_chr_num_str
        #this_chr = this_chr_num_str      #other genome
    this_begin = str(this_window_begin)[-9:]
    this_end = str(this_window_end)[-9:]
    os.system("quickgt " + bamFile + " " + refGenome + " " + this_chr + ":" + this_begin + "-" + this_end + " > tmp.tsv")
    depth_sum = 0
    count = 0
    ins_count = 0
    del_count = 0
    sub_count = 0
    with open('tmp.tsv') as f:
        tsvreader = csv.reader(f, delimiter = '\t')
        for line in tsvreader:
            if line[0] == 'CHR':
                continue
            count = count + 1
            depth_sum = depth_sum + int(line[3])
            if line[9] != "0;0":
                ins_count = ins_count + 1
            elif line[10] != "0;0":
                del_count = del_count + 1
            elif int(line[3]) > (int(line[pileup_tsv_dic[line[2]]].split(";")[0]) + int(line[pileup_tsv_dic[line[2]]].split(";")[1])):
                sub_count = sub_count + 1
    try:
        depth = float(depth_sum) / float(count)
        insRatio = float(ins_count) / float(count)
        delRatio = float(del_count) / float(count)
        subRatio = float(sub_count) / float(count)
    except:
        depth = float(0.1)
        insRatio = float(0.0)
        delRatio = float(0.0)
        subRatio = float(0.0)
    window_SVdepth_arr.append(depth)
    window_SVinsRatio_arr.append(insRatio)
    window_SVdelRatio_arr.append(delRatio)
    window_SVsubRatio_arr.append(subRatio)

snp_WindowsMetaFeature_arr = []
snp_allWindowMetaFeature_file = open('tmp_snp_windowsMetaFeature.csv', 'ab')
for i in range(len(snp_window_arr)):
    if i == len(snp_window_arr) - 1:
        snp_allWindowMetaFeature_file.write(str(snp_window_arr[i][0]) + "," + str(snp_window_arr[i][1]) + "," + str(window_readLen_arr[i]) + "," + str(window_depth_arr[i]) + "," + str(window_sequencingQuality_arr[i]) + "," + str(window_mappingQuality_arr[i]) + "," + str(window_repeatRatio_arr[i]) + "," + str(window_insRatio_arr[i]) + "," + str(window_delRatio_arr[i]) + "," + str(window_subRatio_arr[i]))
    else:
        snp_allWindowMetaFeature_file.write(str(snp_window_arr[i][0]) + "," + str(snp_window_arr[i][1]) + "," + str(window_readLen_arr[i]) + "," + str(window_depth_arr[i]) + "," + str(window_sequencingQuality_arr[i]) + "," + str(window_mappingQuality_arr[i]) + "," + str(window_repeatRatio_arr[i]) + "," + str(window_insRatio_arr[i]) + "," + str(window_delRatio_arr[i]) + "," + str(window_subRatio_arr[i]) + "\n")
    snp_WindowsMetaFeature_arr.append([snp_window_arr[i][0], snp_window_arr[i][1], window_readLen_arr[i], window_depth_arr[i], window_sequencingQuality_arr[i], window_mappingQuality_arr[i], window_repeatRatio_arr[i], window_insRatio_arr[i], window_delRatio_arr[i], window_subRatio_arr[i]])
snp_allWindowMetaFeature_file.close()

sv_WindowsMetaFeature_arr = []
sv_allWindowMetaFeature_file = open('tmp_sv_windowsMetaFeature.csv', 'ab')
for i in range(len(sv_window_arr)):
    if i == len(sv_window_arr) - 1:
        sv_allWindowMetaFeature_file.write(str(sv_window_arr[i][0]) + "," + str(sv_window_arr[i][1]) + "," + str(window_SVreadLen_arr[i]) + "," + str(window_SVdepth_arr[i]) + "," + str(window_length_arr[i]) + "," + str(window_insertSize_arr[i]) + "," + str(window_SVrepeatRatio_arr[i]) + "," + str(window_SVinsRatio_arr[i]) + "," + str(window_SVdelRatio_arr[i]) + "," + str(window_SVsubRatio_arr[i]))
    else:
        sv_allWindowMetaFeature_file.write(str(sv_window_arr[i][0]) + "," + str(sv_window_arr[i][1]) + "," + str(window_SVreadLen_arr[i]) + "," + str(window_SVdepth_arr[i]) + "," + str(window_length_arr[i]) + "," + str(window_insertSize_arr[i]) + "," + str(window_SVrepeatRatio_arr[i]) + "," + str(window_SVinsRatio_arr[i]) + "," + str(window_SVdelRatio_arr[i]) + "," + str(window_SVsubRatio_arr[i]) + "\n")
    sv_WindowsMetaFeature_arr.append([sv_window_arr[i][0], sv_window_arr[i][1], window_SVreadLen_arr[i], window_SVdepth_arr[i], window_length_arr[i], window_insertSize_arr[i], window_SVrepeatRatio_arr[i], window_SVinsRatio_arr[i], window_SVdelRatio_arr[i], window_SVsubRatio_arr[i]])
sv_allWindowMetaFeature_file.close()

os.system("rm tmp*")        #清除临时文件


snp_WindowsMetaFeature_dic = {}
for one_snp_WindowsMetaFeature_arr_ele in snp_WindowsMetaFeature_arr:
    this_metaArr = one_snp_WindowsMetaFeature_arr_ele[2:]
    this_predicted_caller = snp_model.predict([this_metaArr])[0]
    snp_WindowsMetaFeature_dic[str(one_snp_WindowsMetaFeature_arr_ele[0]) + "_" + str(one_snp_WindowsMetaFeature_arr_ele[1])] = this_predicted_caller

sv_WindowsMetaFeature_dic = {}
for one_sv_WindowsMetaFeature_arr_ele in sv_WindowsMetaFeature_arr:
    this_metaArr = one_sv_WindowsMetaFeature_arr_ele[2:]
    this_predicted_caller = sv_model.predict([this_metaArr])[0]
    sv_WindowsMetaFeature_dic[str(one_sv_WindowsMetaFeature_arr_ele[0]) + "_" + str(one_sv_WindowsMetaFeature_arr_ele[1])] = this_predicted_caller


def readIn(iniFileDir):
    iniObject = open(iniFileDir)
    try:
        iniFile = iniObject.read( )
    finally:
        iniObject.close( )
    iniFile_rows = iniFile.split('\n')
    iniFile_rows.pop()
    return iniFile_rows
def get_all_values(dict1, dict2):
  values = set()
  values.update(dict1.values())
  values.update(dict2.values())
  return list(values)
all_values = get_all_values(snp_WindowsMetaFeature_dic, sv_WindowsMetaFeature_dic)
ini_info_dic = {"bcftools": [], "Delly": [], "freebayes": [], "gatk": [], "Manta": [], "Pindel": [], "SvABA": []}
TMBstable_vcf = open(outFileDir, 'ab')
TMBstable_vcf.write("##fileformat=VCFv4.2\n")
for one_all_values_ele in all_values:
    this_dir = ini_info_dir + "/" + one_all_values_ele + ".vcf"
    this_dir_arr = readIn(this_dir)
    ini_info_dic[one_all_values_ele] = this_dir_arr
    for one_this_dir_arr_ele in this_dir_arr:
        if one_this_dir_arr_ele.startswith("##INFO") or one_this_dir_arr_ele.startswith("##FORMAT"):
            TMBstable_vcf.write(one_this_dir_arr_ele + "\n")
        if not one_this_dir_arr_ele.startswith("#"):
            break
TMBstable_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown\n")
for one_snp_WindowsMetaFeature_dic_ele in snp_WindowsMetaFeature_dic:
    this_info = snp_WindowsMetaFeature_dic[one_snp_WindowsMetaFeature_dic_ele]
    this_start = int(one_snp_WindowsMetaFeature_dic_ele.split("_")[0])
    this_end = int(one_snp_WindowsMetaFeature_dic_ele.split("_")[1])
    for one_this_info_dic_ele in ini_info_dic[this_info]:
        if not one_this_info_dic_ele.startswith("#"):
          try:
            this_pos_arr = one_this_info_dic_ele.split("\t")
            this_pos = chromosomesCodingDic[this_pos_arr[0]] + int(this_pos_arr[1])
            if this_pos >= this_start and this_pos <= this_end:
                TMBstable_vcf.write(one_this_info_dic_ele + "\n")
          except:
            continue
for one_sv_WindowsMetaFeature_dic_ele in sv_WindowsMetaFeature_dic:
    this_info = sv_WindowsMetaFeature_dic[one_sv_WindowsMetaFeature_dic_ele]
    this_start = int(one_sv_WindowsMetaFeature_dic_ele.split("_")[0])
    this_end = int(one_sv_WindowsMetaFeature_dic_ele.split("_")[1])
    for one_this_info_dic_ele in ini_info_dic[this_info]:
        if not one_this_info_dic_ele.startswith("#"):
          try:
            this_pos_arr = one_this_info_dic_ele.split("\t")
            this_pos = chromosomesCodingDic[this_pos_arr[0]] + int(this_pos_arr[1])
            if this_pos >= this_start and this_pos <= this_end:
                TMBstable_vcf.write(one_this_info_dic_ele + "\n")
          except:
            continue

