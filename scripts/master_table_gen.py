from genome_comparison import *
from combined import *
import pandas as pd
import os
import csv
import copy
import subprocess
import argparse

genomes = []
PWD = 'test_master_table/test_comparisons'
CHEBI_PATH = 'src/chebi.obo'
SUBSTRATE_TABLE_PATH = 'src/substrateTypes.tsv'
OUTPUT_DIR = 'src/tables'

superfamilies = ['1.A.17', '2.A.7', '1.E.11', '1.E.36', '1.C.12', '1.A.72', '2.A.29', '2.A.66','3.A.3', '1.C.41', '2.A.6']


# generates the columns needed for our Master Dataframe
def master_df_column_gen(directory):
    columns = ['#TCID', 'Acc', 'CE', 'Role', 'hit_tms_no']
    per_genome_columns = ['query', 'q_tms', 'evalue', 'pident', 'qcov', 'scov']
    num_genomes = 0
    items = os.listdir(directory)
    
    global genomes
    genomes = [item for item in items if (os.path.isdir(os.path.join(directory, item)) and 'GCF' in item)]
    
    columns.extend(genomes)
    for genome in genomes:
        columns.extend(per_genome_columns)

    return columns

# declare master_dict as global
master_dict = {}

def ce_role_dict_generation():    
    substrate_dict = get_substrate_data('https://tcdb.org/cgi-bin/substrates/getSubstrates.py')

    input = list(master_dict.keys())

    cleanChebiIDs = get_chebi_id(substrate_dict, input)
    
    role_classes = set(['CHEBI:23888','CHEBI:33281','CHEBI:26672','CHEBI:31432','CHEBI:33229','CHEBI:23357','CHEBI:25212',
                        'CHEBI:23924', 'CHEBI:27780', 'CHEBI:35703', 'CHEBI:37958', 'CHEBI:38161', 'CHEBI:71338', 'CHEBI:62488',
                        'CHEBI:33280', 'CHEBI:25728'])

    graph, idToName, nameToId, idToParent, altIDToID, idToRelationship = oboParse(CHEBI_PATH)
    predecessorDict = terminalPredecessorParse(SUBSTRATE_TABLE_PATH, altIDToID)

    myGroup = {}
    for key in cleanChebiIDs:
        if cleanChebiIDs[key] == 'none':
            myGroup[key] = {}
            myGroup[key]['CE'] = ['none']
            myGroup[key]['Role'] = ['none']
        else:
            myGroup[key] = {}
            myGroup[key]['CE'] = []
            myGroup[key]['Role'] = []
            for id in cleanChebiIDs[key]:
                id = findPrimary(id, altIDToID)#id, id_to_relationship, idToName, altIDToID, graph, idToParent
                idPredecessors = findPredecessor(id,predecessorDict,graph, altIDToID, idToName)
                idRoles = findPredecessor(id,predecessorDict,graph, altIDToID, idToName)
                idRoles = findRole(id,idToRelationship,idToName, altIDToID, graph, idToParent,classes=role_classes)
                if len(idPredecessors) == 0:
                    myGroup[key]['CE'].append(getSubstrateName(id, idToName) + '(' + id + ')')
                else:
                    myGroup[key]['CE'].append(str(idPredecessors) + '-' + getSubstrateName(id, idToName) + '(' + id + ')')
                if len(idRoles) == 0:
                    myGroup[key]['Role'].append(getSubstrateName(id, idToName) + '(' + id + ')')
                else:
                    myGroup[key]['Role'].append(str(idRoles) + '-' + getSubstrateName(id, idToName) + '(' + id + ')')
    
    return myGroup

def master_dict_generation(green_dir):
    getSmithWaterman(PWD)
    parse_sw(PWD)


    
    for filename in os.listdir(green_dir):
        if os.path.isfile(os.path.join(green_dir, filename)):

            file_parts = filename.split('_')[:2]
            genome = '_'.join(file_parts)

            green_df = pd.read_csv(os.path.join(green_dir, filename), delimiter='\t')
            
            for index, row in green_df.iterrows():
                tcid_acc = row['Hit_tcid'] + '-' + row['Hit_xid']

                info_df = None
                try:
                    info_df = getInfoAsRow(genome, tcid_acc)
                except:
                    fix_cmd = f"grep {row['Hit_xid']} ~/db/blastdb/tcdb.faa"
                    correct_tcdb = subprocess.run([fix_cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True).stdout.strip().split('>')[1]
                    tcid_acc = correct_tcdb
                

                if tcid_acc in master_dict and master_dict[tcid_acc]:
                    if genome in master_dict[tcid_acc] and master_dict[tcid_acc][genome]:
                        master_dict[tcid_acc][genome][row['#Query_id']] = {}
                    else:
                        master_dict[tcid_acc][genome] = {row['#Query_id']:{}}
                else: 
                    master_dict[tcid_acc] = {genome: {row['#Query_id']:{}}}
                if info_df is None:
                    master_dict[tcid_acc][genome][row['#Query_id']] = {
                        'CE' : 'N/A',
                        'Role' : 'N/A',
                        'hit_tms_no' : row['Hit_n_TMS'],
                        'query': 'none',
                        'q_tms': row['Query_n_TMS'],
                        'evalue' : 'none',
                        'pident': 'none',
                        'qcov': 'none',
                        'scov': 'none',
                    } 
                else:

                    def which_row(df, query_to_find):
                        for i, r in df.iterrows():
                            if r['hit_id'] == query_to_find:
                                return i
                        
                        return 0

                    idx = which_row(info_df, row['#Query_id'])
                    info_row = info_df.iloc[idx]
               
                    qcov_val = round(float(((info_row['hit_end'] - info_row['hit_start'] + 1) / info_row['hit_len']) * 100), 1)
                    scov_val = round(float(((info_row['query_end'] - info_row['query_start'] + 1) / info_row['query_len']) * 100), 1)


                    master_dict[tcid_acc][genome][row['#Query_id']] = {
                        'CE' : 'N/A',
                        'Role' : 'N/A',
                        'hit_tms_no' : row['Hit_n_TMS'],
                        'query': info_row['hit_id'],
                        'q_tms': row['Query_n_TMS'],
                        'evalue' :info_row['eval'],
                        'pident':info_row['pident'],
                        'qcov': qcov_val,
                        'scov': scov_val,
                    }
    
    ce_role_dict = ce_role_dict_generation()

    for tc_acc in master_dict:
        tcid = tc_acc.split('-')[0]
        ce = ce_role_dict[tcid]['CE']
        role = ce_role_dict[tcid]['Role']
        for genome in master_dict[tc_acc]:
            for q_id in master_dict[tc_acc][genome]:
                master_dict[tc_acc][genome][q_id]['CE'] = ce
                master_dict[tc_acc][genome][q_id]['Role'] = role

NUM_COL_PRE_GENOME = 5
NUM_COL_PER_GENOME = 6

def tsv_generation():
    df_columns = master_df_column_gen(PWD)
    master_dict_generation(PWD + '/greens')
    master_df = pd.DataFrame(columns=df_columns)

    master_array = []
    master_array.append(df_columns)

    for tcid_acc_key in master_dict:
        output_row = [None] * len(df_columns)
        output_row[0] = tcid_acc_key.split('-')[0].strip()
        output_row[1] = tcid_acc_key.split('-')[1].strip()
        

        to_write = False
        all_negative = True

        for i in range(len(genomes)):
            if genomes[i] in master_dict[tcid_acc_key]:

                protein = next(iter(master_dict[tcid_acc_key][genomes[i]]))

                if master_dict[tcid_acc_key][genomes[i]][protein]['q_tms'] == 'none':
                    output_row[NUM_COL_PRE_GENOME+i] = '-'
                else:
                    output_row[NUM_COL_PRE_GENOME+i] = '+'
                    all_negative = False

                min_evalue = float('inf')
                selected_protein = None
                for protein in master_dict[tcid_acc_key][genomes[i]]:
                    protein_vals = master_dict[tcid_acc_key][genomes[i]][protein]
                    evalue = protein_vals.get('evalue', 'none')
                
                    if evalue == 'none':
                        continue

                    evalue = float(evalue)
                    if evalue < min_evalue:
                        min_evalue = evalue
                        selected_protein = protein
                
                if selected_protein != None:
                    to_write = True
                
                
                    first_protein_vals = master_dict[tcid_acc_key][genomes[i]][selected_protein]

                    start_index = NUM_COL_PRE_GENOME+len(genomes)+NUM_COL_PER_GENOME*i-1
                    output_row[start_index+1] = str(first_protein_vals['query']).strip()
                    output_row[start_index+2] = str(first_protein_vals['q_tms']).strip()
                    output_row[start_index+3] = str(first_protein_vals['evalue']).strip()
                    output_row[start_index+4] = str(first_protein_vals['pident']).strip()
                    output_row[start_index+5] = str(first_protein_vals['qcov']).strip()
                    output_row[start_index+6] = str(first_protein_vals['scov']).strip()
                    output_row[2] = ", ".join(first_protein_vals['CE']).replace("'", "").strip()
                    output_row[3] = ", ".join(first_protein_vals['Role']).replace("'", "").strip()
                    output_row[4] = str(first_protein_vals['hit_tms_no']).strip()
            else:
                output_row[NUM_COL_PRE_GENOME+i] = '-'
                start_index = NUM_COL_PRE_GENOME+len(genomes)+NUM_COL_PER_GENOME*i
                for j in range(6):
                    output_row[start_index+j] = 'none'
        if to_write or all_negative:
            master_array.append(output_row)

    with open(f'{OUTPUT_DIR}/master_table.tsv', 'w', newline='') as tsv_output:
        sorted_master_array = sorted(master_array, key=lambda x: str(x[0]))
        writer = csv.writer(tsv_output, delimiter='\t')
        writer.writerows(sorted_master_array)
    
    

def additional_row(prev_row, protein_vals):
        additional_row = prev_row
        additional_row[2] = protein_vals['CE']
        additional_row[3] = protein_vals['Role']
        additional_row[4] = protein_vals['hit_tms_no']
        additional_row[5] = protein_vals['query']
        additional_row[6] = protein_vals['q_tms']
        additional_row[7] = protein_vals['evalue']
        additional_row[8] = protein_vals['pident']
        additional_row[9] = protein_vals['qcov']
        additional_row[10] = protein_vals['scov']
        return additional_row

def genome_tsv_generation():
    master_array = []
    genome_columns = ['#TCID','Acc','CE','Role','hit_tms_no','query','q_tms','evalue','pident','qcov','scov']
    master_array.append(genome_columns)

    for genome in genomes:

        genome_master_array = []

        for tcid_acc_key in master_dict:
            if genome in master_dict[tcid_acc_key]:
                output_row = [None] * len(genome_columns)
                output_row[0] = tcid_acc_key.split('-')[0]
                output_row[1] = tcid_acc_key.split('-')[1]


                temp_arr = []
                for protein in master_dict[tcid_acc_key][genome]:
                    temp = copy.copy(output_row)
                    genome_master_array.append(additional_row(temp, master_dict[tcid_acc_key][genome][protein]))
        
        filename = f"{PWD}/{genome}/test_genome_{genome}.tsv"
        with open(filename, 'w', newline='') as tsv_output:
            writer = csv.writer(tsv_output, delimiter='\t')
            writer.writerows(genome_master_array)
        print('Finished with ' + genome)


def gen_family_tsv(directory, greens_dir):
    family_data = {}
    master_fam_dict = {}
    items = os.listdir(directory)
    genome_content = {}
    total_families = []
    genomes = [item for item in items if (os.path.isdir(os.path.join(directory, item)) and 'GCF' in item)]

    fam_df = pd.DataFrame(columns=genomes)
    for filename in os.listdir(greens_dir):
        file_parts = filename.split('_')[:2]
        genome = '_'.join(file_parts)

        green_df = pd.read_csv(os.path.join(greens_dir, filename), delimiter='\t')
        
        genome_content[genome] = {}
        for index, row in green_df.iterrows():
            tcid = row['Hit_tcid']

            temp = tcid.split('.')
            # define a family either on the first 3 characters of tcid or first 4
            fam = temp[0] + '.' + temp[1] + '.' + temp[2]
            if fam in superfamilies:
                fam = temp[0] + '.' + temp[1] + '.' + temp[2] + '.' + temp[3] 

            if fam not in total_families:
                total_families.append(fam)
            
            if fam not in genome_content[genome]:
                genome_content[genome][fam] = 1
            else:
                genome_content[genome][fam] += 1

    for g in genome_content:
        family_data[g] = []

        for f in total_families:
            if f in genome_content[g]:
                family_data[g].append(1)
            else:
                family_data[g].append(0)

    fam_percents = {}
    for f in total_families:
        fam_percents[f] = []
        
        for g in genome_content:
            if f in genome_content[g]:
                fam_percents[f].append(genome_content[g][f])
            else:
                fam_percents[f].append(0)
    


    for family in fam_percents:
        max_val = max(fam_percents[family])
        
        new_list = []
        for val in fam_percents[family]:
            percent = round(float(val/max_val), 2)
            new_list.append(percent)
        
        fam_percents[family] = new_list

    


    percent_df = pd.DataFrame.from_dict(fam_percents, orient='index')
    percent_df.columns = genomes
    percent_df.to_csv(f'{OUTPUT_DIR}/family_percent.tsv')
    
    fam_df = pd.DataFrame(family_data)
    fam_df.index = total_families
    fam_df.to_csv(f'{OUTPUT_DIR}/family_matrix.tsv')

def main():
    parser = argparse.ArgumentParser(description='Master Table Generation Script')
    parser.add_argument('-p', '--PWD', type=str, help='Greens Information Directory')
    parser.add_argument('-c', '--chebi', type=str, help='Chebi obo file location')
    parser.add_argument('-s', '--substrates', type=str, help='Substrate table file location')
    parser.add_argument('-o', '--outDir', type=str, help='Output directory')

    args = parser.parse_args()

    global PWD
    global CHEBI_PATH
    global SUBSTRATE_TABLE_PATH
    global OUTPUT_DIR

    if args.PWD is not None:
        PWD = args.PWD
    if args.chebi is not None:
        CHEBI_PATH = args.chebi
    if args.substrates is not None:
        SUBSTRATE_TABLE_PATH = args.substrates
    if args.outDir is not None:
        OUTPUT_DIR = args.outDir
    
    master_dict_generation(PWD + '/greens')
    master_df_column_gen(PWD)
    #genome_tsv_generation()
            

    tsv_generation()
    #gen_family_tsv(PWD, PWD + '/greens')

main()

