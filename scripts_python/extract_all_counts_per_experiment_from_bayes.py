import numpy as np
import pandas as pd
import re

#Functions

def remove_underscore_column_names(dataset):
    dataset_renamed = dataset.rename(columns=lambda x: re.sub(r"_.*", "",x))
    return dataset_renamed

def generate_samplefile(samplefilename):
    samplefile =[]
    with open(samplefilename) as infile:
        for line in infile: #remove everything after the first underscoe and remove newline
            samplefile.append(re.sub(r"_.*", "", str(line.rstrip("\n"))))
    return samplefile

def save_df_for_samplefile(samplefile, samplefilename):
    df_samplefile = df_bayes[samplefile]
    ###df_samplefile = remove_3quarters_NAs(df_samplefile)
    df_samplefile = df_samplefile.fillna(0)
    df_samplefile.to_csv(str(samplefilename)+'_counts.csv')

def remove_3quarters_NAs(dataframe):
    df_with_NAs = dataframe.replace('0', np.nan)
    thresh = int(len(df_with_NAs.columns)/4*3)
    NA_removed =df_with_NAs.dropna(axis='index', thresh=thresh)
    return NA_removed


#script

df_bayes = pd.DataFrame.from_csv('/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/bayesian_correction_10Nov/bayesian_corrected_new.csv', header=0, sep=',', index_col=0)
df_bayes = remove_underscore_column_names(df_bayes)

samplefilenameBatch3 = '/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_Batch3.txt'
samplefilenameBatch4 = '/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_Batch4.txt'
samplefilenameBatch5 = '/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_Batch5.txt'
samplefilenameBatch8 = '/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_Batch8.txt'
samplefilenameSS0209 = '/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_SS0209.txt'
samplefilenameSS2608 = '/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_SS2608.txt'
samplefilename_allp300 ='/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_allp300.txt'
samplefilename_allSet7 ='/Users/anna_koeferle/Documents/UCL/HiSeqRun_Sept2015/gRNA_counts_bayesian_correction/Samplefiles/Samplefile_allSet7.txt'

samplefileBatch3 = generate_samplefile(samplefilenameBatch3)
samplefileBatch4 = generate_samplefile(samplefilenameBatch4)
samplefileBatch5 = generate_samplefile(samplefilenameBatch5)
samplefileBatch8 = generate_samplefile(samplefilenameBatch8)
samplefileSS0209 = generate_samplefile(samplefilenameSS0209)
samplefileSS2608 = generate_samplefile(samplefilenameSS2608)
allp300samples = generate_samplefile(samplefilename_allp300)
allSet7samples = generate_samplefile(samplefilename_allSet7)

#separate set7 and p300s if together in same sample
samplefileSS2608_p300 = samplefileSS2608[0:4]
samplefileSS2608_Set7 = samplefileSS2608[0:2] + samplefileSS2608[4:6]
sampleSS0209_p300 = samplefileSS0209[4:8] + samplefileSS0209[10:12]
sampleSS0209_Set7 = samplefileSS0209[0:1] + samplefileSS0209[9:12] + samplefileSS0209[13:15]
sampleBatch3_Set7 = samplefileBatch3[0:3] + samplefileBatch3[6:8]

save_df_for_samplefile(samplefileBatch4, 'Batch4_p300')
save_df_for_samplefile(samplefileBatch5, 'Batch5_Set7')
save_df_for_samplefile(samplefileBatch8, 'Batch8_p300')
save_df_for_samplefile(samplefileSS2608_p300, 'BatchSS2608_p300')
save_df_for_samplefile(samplefileSS2608_Set7, 'BatchSS2608_Set7')
save_df_for_samplefile(sampleSS0209_p300, 'BatchSS0209_p300')
save_df_for_samplefile(sampleSS0209_Set7, 'BatchSS0209_Set7')
save_df_for_samplefile(sampleBatch3_Set7, 'Batch3_Set7')
save_df_for_samplefile(allp300samples, 'All_p300')
save_df_for_samplefile(allSet7samples, 'All_Set7')

print "DONE!"
