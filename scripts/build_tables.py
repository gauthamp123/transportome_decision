import numpy as np
import pandas as pd
import csv
import os

list_of_families = ['1.A.134', '1.W.5', '9.A.32', '9.B.70', '9.B.111', '9.B.325', '9.B.354', '9.B.402', '1.E.29', '1.W.9', '2.A.7.18', '3.A.3.23', '3.A.13', '9.B.323', '1.E.6', '1.W.1', '9.B.270', '9.B.302', '1.C.3', '1.E.11.1', '1.E.31', '2.A.7.5', '9.B.446', '1.A.2', '1.B.89', '1.B.95', '1.C.36', '1.A.44', '1.C.23', '1.E.5', '1.M.12', '1.M.13', '9.B.69', '1.B.59', '1.C.7', '2.A.7.9', '9.B.50', '1.B.16', '2.A.29.14', '2.A.66.10', '8.A.128', '9.B.155', '9.B.171', '1.B.162', '2.A.29.11', '1.C.20', '1.C.37', '2.A.66.5', '2.A.135', '1.E.18', '9.B.412', '1.A.132', '1.B.4', '1.B.79', '1.B.81', '1.C.65', '2.A.7.30', '2.A.12', '5.B.9', '8.A.82', '8.A.232', '9.B.30', '9.B.40', '9.B.110', '9.B.117', '9.B.353', '9.B.358']
superfamilies = ['1.A.17', '2.A.7', '1.E.11', '1.E.36', '1.C.12', '1.A.72', '2.A.29', '2.A.66','3.A.3', '1.C.41', '8.A.180', '2.A.6']

list_of_families_neg = ['2.A.73', '8.A.47', '1.D.194', '9.B.122', '4.C.2', '1.C.12.1', '2.A.130', '2.A.131', '1.B.77', '9.B.176', '9.B.392', '2.A.116', '9.A.10', '1.A.72.3', '1.A.72.4', '1.A.72.5', '1.E.48', '1.M.7', '2.A.7.22', '1.C.4', '1.C.96', '5.B.15', '1.B.62', '1.B.85', '1.C.131', '2.A.7.36', '9.B.23', '9.B.96', '9.A.31', '9.B.139', '9.B.118', '9.B.164']

master_table = pd.read_csv('new_sorted_master_table.tsv', sep='\t')

def starts_with_family_prefix(value):
    for family in list_of_families:
        if value.startswith(family):
            return True
    return False

master_table_covid = master_table[master_table['#TCID'].apply(starts_with_family_prefix)]


def gen_negative_control_tsv(directory, greens_dir):
    family_data = {}
    master_fam_dict = {}
    items = os.listdir(directory)
    genome_content = {}
    total_proteins = []
    genomes = [item for item in items if (os.path.isdir(os.path.join(directory, item)) and 'GCF' in item)]
    substrate_info = {}

    covid_df = pd.DataFrame(columns=genomes)
    for filename in os.listdir(greens_dir):
        file_parts = filename.split('_')[:2]
        genome = '_'.join(file_parts)
 
        green_df = pd.read_csv(os.path.join(greens_dir, filename), delimiter='\t')
        
        genome_content[genome] = {}
        for index, row in green_df.iterrows():
            tcid = row['Hit_tcid']
            acc = row['Hit_xid']
            query = row['#Query_id']
            tcid_acc = tcid + '-' + acc


            temp = tcid.split('.')
            # define a family either on the first 3 characters of tcid or first 4
            fam = temp[0] + '.' + temp[1] + '.' + temp[2]
            if fam in superfamilies:
                fam = temp[0] + '.' + temp[1] + '.' + temp[2] + '.' + temp[3] 
 
            if fam in list_of_families_neg:
                genome_content[genome][tcid_acc] = query
                total_proteins.append(tcid_acc)
                substrate_info[tcid_acc] = row['Predicted_Substrate']
            
            
    for g in genome_content:
        family_data[g] = []

        for f in total_proteins:
            if f in genome_content[g]:
                family_data[g].append(genome_content[g][f])
            else:
                family_data[g].append('-')

    neg_df = pd.DataFrame(family_data)
    neg_df.index = total_proteins
    values_list = [substrate_info[row] if row in substrate_info else None for row in neg_df.index]

    # Add a new column using the values list
    #covid_df['Substrate'] = values_list
    neg_df.insert(0, 'Substrate', values_list)
    neg_df = neg_df.rename_axis('Protein')
    neg_df = neg_df.drop_duplicates()
    
    neg_df.to_csv('negcontrol_specific_fams.tsv')

#master_table_covid.to_csv('test.tsv')
 
def gen_covid_tsv(directory, greens_dir):
    family_data = {}
    master_fam_dict = {}
    items = os.listdir(directory)
    genome_content = {}
    total_proteins = []
    genomes = [item for item in items if (os.path.isdir(os.path.join(directory, item)) and 'GCF' in item)]
    substrate_info = {}

    covid_df = pd.DataFrame(columns=genomes)
    for filename in os.listdir(greens_dir):
        file_parts = filename.split('_')[:2]
        genome = '_'.join(file_parts)
 
        green_df = pd.read_csv(os.path.join(greens_dir, filename), delimiter='\t')
        
        genome_content[genome] = {}
        for index, row in green_df.iterrows():
            tcid = row['Hit_tcid']
            acc = row['Hit_xid']
            query = row['#Query_id']
            tcid_acc = tcid + '-' + acc


            temp = tcid.split('.')
            # define a family either on the first 3 characters of tcid or first 4
            fam = temp[0] + '.' + temp[1] + '.' + temp[2]
            if fam in superfamilies:
                fam = temp[0] + '.' + temp[1] + '.' + temp[2] + '.' + temp[3] 
 
            if fam in list_of_families:
                genome_content[genome][tcid_acc] = query
                total_proteins.append(tcid_acc)
                substrate_info[tcid_acc] = row['Predicted_Substrate']
            
            
    for g in genome_content:
        family_data[g] = []

        for f in total_proteins:
            if f in genome_content[g]:
                family_data[g].append(genome_content[g][f])
            else:
                family_data[g].append('-')

    covid_df = pd.DataFrame(family_data)
    print(covid_df.columns)
    covid_df.index = total_proteins
    values_list = [substrate_info[row] if row in substrate_info else None for row in covid_df.index]

    # Add a new column using the values list
    #covid_df['Substrate'] = values_list
    covid_df.insert(0, 'Substrate', values_list)
    covid_df = covid_df.rename_axis('Protein')
    covid_df = covid_df.drop_duplicates()
    
    covid_df.to_csv('covid_specific_proteins.tsv')




    
gen_covid_tsv('test_comparisons', 'test_comparisons/greens')
gen_negative_control_tsv('test_comparisons', 'test_comparisons/greens')