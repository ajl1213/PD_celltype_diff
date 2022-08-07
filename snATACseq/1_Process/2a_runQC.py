#!/home/ajl1213/anaconda2/bin/python


import os



DIR=os.getcwd()

sampleList=[

]


os.system('mkdir '+DIR+'/Plots')


for sampleID in sampleList:
    print sampleID
    runLine='/home/ajl1213/anaconda2/envs/r_env/bin/R CMD BATCH --no-save --no-restore \'--args '+sampleID+'\' '+DIR+'/2b_sampleQC.R '+DIR+'/Plots/'+sampleID+'.sampleQC.Rout'
    print runLine
    os.system(runLine)






