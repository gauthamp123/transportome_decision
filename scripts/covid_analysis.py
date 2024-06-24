import numpy as np
import pandas as pd
import os
import csv

fam_matrix = pd.read_csv('family_matrix.tsv')

covid_genomes = ['GCF_001865765.1', 'GCF_000769635.1', 'GCF_006094375.1', 'GCF_000013425.1', 'GCF_017310215.1', 'GCF_000959505.1', 'GCF_009648975.1', 'GCF_900660445.2', 'GCF_003006415.1', 'GCF_900637655.1', 'GCF_000961215.1', 'GCF_016576965.1', 'GCF_003287895.1', 'GCF_002222595.2', 'GCF_015159595.1', 'GCF_010669205.1', 'GCF_001689125.2']
neg_control = ['GCF_000144405.1','GCF_000240185.1','GCF_001558775.1','GCF_002055515.1','GCF_003367905.1','GCF_008632635.1','GCF_008693705.1','GCF_016698705.1','GCF_016728665.1','GCF_020405345.1','GCF_023520795.1','GCF_900660465.1']

exclusive_to_covid = []
print('In COVID but not in negative control')
for index, row in fam_matrix.iterrows():
    # check all covid genomes
    in_covid = False
    for genome in covid_genomes:
        if row[genome] == 1:
            in_covid = True

    in_control = False
    for genome in neg_control:
        if row[genome] == 1:
            in_control = True
    
    if in_covid == True and in_control == False:
        exclusive_to_covid.append(row['family'])


print(exclusive_to_covid)
print('---------------')

print('In negative control but absent in COVID')
exclusive_to_negative_control = []
for index, row in fam_matrix.iterrows():
    # check all covid genomes
    in_neg_control = False
    for genome in neg_control:
        if row[genome] == 1:
            in_neg_control = True
    
    in_covid_genome = False
    if not in_neg_control:
        continue
    else:
        for genome in covid_genomes:
            if row[genome] == 1:
                in_covid_genome = True
    
    if in_covid_genome == False and in_neg_control == True:
        exclusive_to_negative_control.append(row['family'])



print(exclusive_to_negative_control)
