import numpy as np
import xlmhg
import pandas as pd
import os
import csv
from numpy import genfromtxt
from plotly.offline import plot



dirPath = "/Bigdata/Dropbox (Technion Dropbox)/Rina_Benel/Home/Rina/AS_linc/cardiac_diff/results/salmon10.1_results/Correlation/lincRNA_lincRNAnetwork/XLmhg/SingleModuleUpdatedAllFunc"

for filename in os.listdir(dirPath):
    if filename.endswith('.csv'):
        with open(os.path.join(dirPath, filename)) as f:
            content = genfromtxt(f, delimiter=",", dtype=None, skip_header=1)
            MyV = np.uint(content)
            #print(MyV)


            #v = np.uint([1,0,1,1,0,1] + [0]*12 + [1,0])
            #X = 3
            #L = 10

            stat, cutoff, pval = xlmhg.xlmhg_test(MyV)

            print(filename)
            print('Test statistic: %.3f' % stat)
            print('Cutoff: %d' % cutoff)
            print('P-value: %.3f' % pval)

            N = MyV.size
            indices = np.uint16(np.nonzero(MyV)[0])

            result = xlmhg.get_xlmhg_test_result(N, indices)

            fig = xlmhg.get_result_figure(result, show_inset=False)

            plot(fig, filename= os.path.join(dirPath, filename) + 'PlotParamaters' + 'plotFigure.html')

