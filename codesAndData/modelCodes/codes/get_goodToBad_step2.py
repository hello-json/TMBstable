# -*- coding: UTF-8 -*
import argparse
import os

'''
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--m_shell_dir', help='the mulation_shell_dir', required=True, type=str)
parser.add_argument('-c', '--call_file_dirs', help='the calling_file_dirs, 用分号分割', required=True, type=str)
    
args = parser.parse_args()
answerFileDir = args.m_shell_dir
callFileDirs = args.call_file_dirs
'''
#answerFileDir = "/mnt/GSDcreator/A_batchMuShells/20230701181614_290_100_264_0800_2410_190000.sh"
#callFileDirs = "/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/bcftools/bcftools_result/20230701181614/20230701181614_bcftools.filter.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Delly/Delly_result/20230701181614/20230701181614_Delly.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/freebayes/freebayes_result/20230701181614/20230701181614_freebayes.filter5.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/gatk/gatk_result/20230701181614/20230701181614_gatk.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Manta/Manta_result/20230701181614/results/variants/candidateSmallIndels.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Manta/Manta_result/20230701181614/results/variants/candidateSV.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Manta/Manta_result/20230701181614/results/variants/diploidSV.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Pindel/Pindel_result/20230701181614/20230701181614_D.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Pindel/Pindel_result/20230701181614/20230701181614_INV.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Pindel/Pindel_result/20230701181614/20230701181614_LI.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Pindel/Pindel_result/20230701181614/20230701181614_SI.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Pindel/Pindel_result/20230701181614/20230701181614_TD.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/SvABA/SvABA_result/20230701181614/20230701181614.svaba.indel.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/SvABA/SvABA_result/20230701181614/20230701181614.svaba.sv.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/SvABA/SvABA_result/20230701181614/20230701181614.svaba.unfiltered.indel.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/SvABA/SvABA_result/20230701181614/20230701181614.svaba.unfiltered.sv.vcf"
        
chromosomesCodingDic = {"chr1": 1000000000, "chr2": 2000000000, "chr3": 3000000000, "chr4": 4000000000, "chr5": 5000000000, "chr6": 6000000000, "chr7": 7000000000, "chr8": 8000000000, "chr9": 9000000000, "chr10": 10000000000, "chr11": 11000000000, "chr12": 12000000000, "chr13": 13000000000, "chr14": 14000000000, "chr15": 15000000000, "chr16": 16000000000, "chr17": 17000000000, "chr18": 18000000000, "chr19": 19000000000, "chr20": 20000000000, "chr21": 21000000000, "chr22": 22000000000, "chrX": 23000000000, "chrY": 24000000000, "1": 1000000000, "2": 2000000000, "3": 3000000000, "4": 4000000000, "5": 5000000000, "6": 6000000000, "7": 7000000000, "8": 8000000000, "9": 9000000000, "10": 10000000000, "11": 11000000000, "12": 12000000000, "13": 13000000000, "14": 14000000000, "15": 15000000000, "16": 16000000000, "17": 17000000000, "18": 18000000000, "19": 19000000000, "20": 20000000000, "21": 21000000000, "22": 22000000000, "X": 23000000000, "Y": 24000000000}

def get_evaluation_value(answerFileDir, callFileDirs):
    answerObject = open(answerFileDir)
    try:
         answerFile = answerObject.read( )
    finally:
         answerObject.close( )
    answerFile_rows = answerFile.split('\n')
    answerFile_arr = answerFile_rows[0].split(' -')
    #print answerFile_arr
    snpAnswerArr = answerFile_arr[4].split(' ')[1].split(':')
    delAnswerArr = answerFile_arr[5].split(' ')[1].split(':')
    insAnswerArr = answerFile_arr[6].split(' ')[1].split(':')
    trAnswerArr = answerFile_arr[7].split(' ')[1].split(':')
    invAnswerArr = answerFile_arr[8].split(' ')[1].split(':')
    allAnswerArr = snpAnswerArr + delAnswerArr + insAnswerArr + trAnswerArr + invAnswerArr
    #print allAnswerArr
    answerPosArr = []
    for one_answer in allAnswerArr:
        one_answer_chr = one_answer.split(',')[1]
        one_answer_ini_pos = int(one_answer.split(',')[2])
        if one_answer_chr in chromosomesCodingDic:
            answerPosArr.append(chromosomesCodingDic[one_answer_chr] + one_answer_ini_pos)
    answerPosArr = list(set(answerPosArr))
    answerPosArr.sort()
    #print answerPosArr
    
    callFileDir_arr = callFileDirs.split(";")
    callPosArr = []
    for callFileDir in callFileDir_arr:
        try:
            callObject = open(callFileDir)
        except:
            callPosArr = []
            break
        try:
            callFile = callObject.read( )
        finally:
            callObject.close( )
        callFile_rows = callFile.split('\n')
        callFile_rows.pop()
        #print callFile_rows
        for one_callFile_row in callFile_rows:
            if not one_callFile_row.startswith("#"):
                one_call_chr = one_callFile_row.split('\t')[0]
                one_call_ini_pos = int(one_callFile_row.split('\t')[1])
                if one_call_chr in chromosomesCodingDic:
                    callPosArr.append(chromosomesCodingDic[one_call_chr] + one_call_ini_pos)
    callPosArr = list(set(callPosArr))
    callPosArr.sort()
    #print callPosArr
       
   
    TP = 0
    call_in_ref = []
    for one_callPos in callPosArr:
        for one_answerPos in answerPosArr:
            if abs(one_answerPos - one_callPos) < 3:
                call_in_ref.append(one_answerPos)
                TP = TP + 1
                break
    print TP
   
    TP = 0
    ref_in_call = []
    for one_answerPos in answerPosArr:
        for one_callPos in callPosArr:
            if abs(one_answerPos - one_callPos) < 1:
                ref_in_call.append(one_callPos)
                TP = TP + 1
                break
    if TP == 0:
        for one_answerPos in answerPosArr:
            for one_callPos in callPosArr:
                if abs(one_answerPos - one_callPos) < 10000:
                    ref_in_call.append(one_callPos)
                    TP = TP + 1
                    break
    
     
    call_len = 0
    answerPosFlagArr = [0] * len(answerPosArr)
    for one_callPos in callPosArr:
        flag = 0
        for i in range(len(answerPosArr)):
            if abs(answerPosArr[i] - one_callPos) < 10:
                answerPosFlagArr[i] = 1
                flag = 1
                break
        if flag == 0:
            call_len = call_len + 1
    for one_answerPosFlag in answerPosFlagArr:
        if one_answerPosFlag == 1:
            call_len = call_len + 1
    #print call_len
    #call_len = len(callPosArr)

    answer_len = len(answerPosArr)
    
    try:
      
        return 1.0 - (float(call_len) - float(TP)) / float(call_len), 1.0 - (float(answer_len) - float(TP)) / float(answer_len)
    except:
        return 1.0 - float(1), 1.0 - float(1)

def findAllDir(base):
    for root, ds, fs in os.walk(base):
        #for f in fs:
        for f in ds:
            yield f

def findAllFile(base):
    for root, ds, fs in os.walk(base):
        #for f in fs:
        for f in fs:
            yield f

def get_callFileDirs_arr(the_sample_name_prefix):
    bcftools_fileDirs = "/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/bcftools/bcftools_result/" + the_sample_name_prefix + "/" + the_sample_name_prefix + "_bcftools.filter.vcf"
    Delly_fileDirs = "/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Delly/Delly_result/" + the_sample_name_prefix + "/" + the_sample_name_prefix + "_Delly.vcf"
    freebayes_fileDirs = "/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/freebayes/freebayes_result/" + the_sample_name_prefix + "/" + the_sample_name_prefix + "_freebayes.filter5.vcf"
    gatk_fileDirs = "/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/gatk/gatk_result/" + the_sample_name_prefix + "/" + the_sample_name_prefix + "_gatk.vcf"
    Manta_fileDirs = "/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Manta/Manta_result/" + the_sample_name_prefix + "/results/variants/candidateSmallIndels.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Manta/Manta_result/" + the_sample_name_prefix + "/results/variants/candidateSV.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Manta/Manta_result/" + the_sample_name_prefix + "/results/variants/diploidSV.vcf"
    Pindel_fileDirs = "/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Pindel/Pindel_result/" + the_sample_name_prefix + "/" + the_sample_name_prefix + "_D.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Pindel/Pindel_result/" + the_sample_name_prefix + "/" + the_sample_name_prefix + "_INV.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Pindel/Pindel_result/" + the_sample_name_prefix + "/" + the_sample_name_prefix + "_LI.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Pindel/Pindel_result/" + the_sample_name_prefix + "/" + the_sample_name_prefix + "_SI.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/Pindel/Pindel_result/" + the_sample_name_prefix + "/" + the_sample_name_prefix + "_TD.vcf"
    SvABA_fileDirs = "/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/SvABA/SvABA_result/" + the_sample_name_prefix + "/" + the_sample_name_prefix + ".svaba.indel.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/SvABA/SvABA_result/" + the_sample_name_prefix + "/" + the_sample_name_prefix + ".svaba.sv.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/SvABA/SvABA_result/" + the_sample_name_prefix + "/" + the_sample_name_prefix + ".svaba.unfiltered.indel.vcf;/mnt/doctoralPeriod/DMADforShortReadsPoint/result/mData/SvABA/SvABA_result/" + the_sample_name_prefix + "/" + the_sample_name_prefix + ".svaba.unfiltered.sv.vcf"
    callFileDirs_arr = [bcftools_fileDirs, Delly_fileDirs, freebayes_fileDirs, gatk_fileDirs, Manta_fileDirs, Pindel_fileDirs, SvABA_fileDirs]
    return callFileDirs_arr

def get_caller_goodToBad(the_sample_calling_result_arr):
    caller_goodToBad_arr = []
    call_name_arr = ["bcftools", "Delly", "freebayes", "gatk", "Manta", "Pindel", "SvABA"]
    the_sample_calling_result_dic = dict(zip(call_name_arr, the_sample_calling_result_arr))
    for i in range(1000):
        if len(the_sample_calling_result_dic) != 0:
            ans = max(the_sample_calling_result_dic, key = lambda x: the_sample_calling_result_dic[x])
            caller_goodToBad_arr.append(ans)
            the_sample_calling_result_dic.pop(ans)
        else:
            break
    return caller_goodToBad_arr

caller_FPR_goodToBad_file = open('mData_allCaller_FPR_goodToBad.csv', 'wb')
caller_FNR_goodToBad_file = open('mData_allCaller_FNR_goodToBad.csv', 'wb')
for the_sample_name in findAllDir("/mnt/GSDcreator/A_batchMuSamples/"):
    this_bam_file = "/mnt/GSDcreator/A_batchMuSamples/" + the_sample_name + "/mmm.bam"
    fsize = os.path.getsize(this_bam_file)
    if fsize != 0:
        print the_sample_name
        the_sample_calling_result_FPR_arr = []
        the_sample_calling_result_FNR_arr = []
        for the_sample_shell in findAllFile("/mnt/GSDcreator/A_batchMuShells/"):
            if the_sample_shell.startswith(the_sample_name):
                this_shell_file = "/mnt/GSDcreator/A_batchMuShells/" + the_sample_shell
                break
        callFileDirs_arr = get_callFileDirs_arr(the_sample_name.split("_")[0])
        for one_callFileDirs in callFileDirs_arr: 
            this_FRP, this_FNR = get_evaluation_value(this_shell_file, one_callFileDirs)
            #print this_FRP, this_FNR
            the_sample_calling_result_FPR_arr.append(this_FRP)
            the_sample_calling_result_FNR_arr.append(this_FNR)
        #print sorted(the_sample_calling_result_FPR_arr), sorted(the_sample_calling_result_FNR_arr)
        caller_FPR_goodToBad_arr = get_caller_goodToBad(the_sample_calling_result_FPR_arr)
        caller_FNR_goodToBad_arr = get_caller_goodToBad(the_sample_calling_result_FNR_arr)
        #if caller_FPR_goodToBad_arr[0] == "freebayes" or caller_FNR_goodToBad_arr[0] == "freebayes":
        #    print caller_FPR_goodToBad_arr, caller_FNR_goodToBad_arr
        one_FPR_line_content = the_sample_name
        for one_caller_FPR_goodToBad in caller_FPR_goodToBad_arr:
            one_FPR_line_content = one_FPR_line_content + "," + one_caller_FPR_goodToBad
        caller_FPR_goodToBad_file.write(one_FPR_line_content + "\n")
        one_FNR_line_content = the_sample_name 
        for one_caller_FNR_goodToBad in caller_FNR_goodToBad_arr:
            one_FNR_line_content = one_FPR_line_content + "," + one_caller_FNR_goodToBad
        caller_FNR_goodToBad_file.write(one_FNR_line_content + "\n")
caller_FPR_goodToBad_file.close()
caller_FNR_goodToBad_file.close()

