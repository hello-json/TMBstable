# -*- coding: UTF-8 -*-

import os

def findAllFile(base):
    for root, ds, fs in os.walk(base):
        #for f in fs:
        for f in ds:
            yield f

i = 0
for the_sample_name in findAllFile("/A_batchSimuSamples/"):
    #print the_sample_name
    this_bam_file = "/ A_batchSimuSamples/" + the_sample_name + "/mmm.bam"
    fsize = os.path.getsize(this_bam_file)
    #print fsize
    if fsize != 0:
        i = i + 1
        print i
        os.system("python oneSample_featureExtraction_writeInFileVersion.py -b " + this_bam_file)


