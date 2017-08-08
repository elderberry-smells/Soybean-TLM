# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 13:16:54 2017
Version 1.0 of the Soybean Trait Linked Marker Tool
@author: U590135
"""
import csv
import pandas as pd
from pandas import ExcelWriter
import glob
import os
import time


haplotype_panel = {
               'rps1_haplotypes' : [
                   # the following allele combinations will produce 'Trait'
                   {'DAS_Gm03_4458273_G_A Zygosity Call': 'G:G',
                    'DAS_Gm03_3838302_A_T Zygosity Call': 'A:A'},
                
                   {'DAS_Gm03_4458273_G_A Zygosity Call': 'A:A',
                    'DAS_Gm03_3838302_A_T Zygosity Call': 'T:T'}
                    ],

                'rps3a_haplotypes' : [
                    # the following allele combinations will produce 'Trait'
                    {'Gm13_28913825_G_A Zygosity Call': 'G:G',
                     'DAS_Gm13_29375136_A_G Zygosity Call': 'G:G'},
                     ],

                'rhg1_haplotypes' : [
                   # the following allele combinations will produce 'Trait'
                   {'DAS_Gm18_1607524_A_C Zygosity Call': 'C:C',
                    'DAS_Gm18_1634714_A_G Zygosity Call': 'G:G',
                    'DAS_Gm18_1686081_A_G Zygosity Call': 'G:G'}
                    ],

                'scc_3_haplotype' : [
                   # the following allele combinations will produce 'Trait'
                   {'Gm03_46333142_G_A Zygosity Call': 'G:G',
                    'Gm03_46508111_A_C Zygosity Call': 'A:A'}
                   ],
                 
                 'scc_11_haplotype' : [
                    # the following allele combinations will produce 'Trait'
                    {'Gm11_37179389_G_T Zygosity Call': 'G:G',
                     'NCSB_002654 Zygosity Call': 'A:A'}
                    ],
                    
                # if adding new panels put them below this line similar to example
                    
                 'example: new panel name here' : [
                    # add in the haplotypes that would produce 'Trait below this line
                    {'example: assay 1 name': 'example: allele combination here',
                     'example: assay 2 name': 'example: allele combination here'},
                     
                     # if there are more than 1 haplotype that produces trait, add comma
                     # at end of the 2nd curly bracket (like above) and make another set like below.
                     # repeat as needed with more sets.  put comma at end of square bracket like above.
                     
                    {'example: assay 1 name': 'example: 2nd combination here',
                     'example: assay 2 name': 'example: 2nd combination here'}
                     ]
                }

def indexAssays(panel, headers):
    """read through Kraken file and return an index of where each panel is located in file.
    If there are additions to any panel or a new panel, they must be added here to the 
    haplotype_panel dictionary
    """
    global haplotype_panel
    
    panel_index = []
    index_assay = haplotype_panel[panel][0]

    for key in index_assay:
        # noinspection PyBroadException
        try:
            if headers.index(key) in panel_index:  # if panel has been appended already, skip
                continue
            else:
                panel_index.append(headers.index(key))  # get index of panel (if exists)
        except:
            continue
        
    return panel_index
                    
        
def getCall(dic, panel):
    '''determine the samples call for each panel in kraken file'''
    global haplotype_panel
    dic_len = len(dic)
    panel_len = len(haplotype_panel[panel])
    
        
    seg = 0
    dupe = 0
    no_data = 0
    trait = 0
    
    #determine if the combination of assays is the Trait haplotype
    if panel_len == 1:
        for assay, allele in dic.items():
            if (assay, allele) in haplotype_panel[panel][0].items():
                trait += 1
        
    else:
        if panel_len == 2:
            for assay, allele in dic.items():
                if (assay, allele) in haplotype_panel[panel][0].items():
                    trait += 1
            if trait != dic_len:
                trait = 0
            
            for assay, allele in dic.items():
                if (assay, allele) in haplotype_panel[panel][1].items():
                    trait += 1
                                
                
            # determine if there is a seg/no call or no data in the panel
    for assay, allele in dic.items():
        if allele == 'No Data':
            no_data += 1
        elif allele == 'Empty':
            no_data += 1
        elif allele == 'Unused':
            no_data += 1
        elif allele == 'Dupe':
            dupe += 1
        elif allele[0] != allele[-1]: 
            seg += 1
             
    # using the counts of each type, return the call        
    if no_data > 0:
        return 'No Data'
    elif dupe > 0:
        return 'No Call'
    elif seg > 0:
        return 'Seg'
    elif trait == dic_len:
        return 'Trait'
    else:
        return 'Wildtype'
    


# ------------------------------------- finding the files in the folder ----------------------------
'''reading the kraken csv file with pandas and converting the data set into individual dataframes'''

path = r'C:\Users\U590135\.spyder-py3\projects\Soybean TLM'
xlextension = 'xlsx'
extension = 'csv'
os.chdir(path)
xlsx_result = [xl for xl in glob.glob('*.{}'.format(xlextension))]

# convert the xlsx files into csv files
xl_files = []
for xlsx_files in xlsx_result:
    if '_report' in xlsx_files:
        continue
    else:
        xl_files.append(xlsx_files)

for i in xl_files:
    df_temp = pd.read_excel(i, index=False, encoding='utf-8')
    get_name = i.find('.')
    csv_name = i[:get_name] + '.csv'
    df_temp.to_csv(csv_name, index=False, encoding='utf-8')
    os.remove(i)  # delete the original XLSX file so there are no duplicates

#  get a list of all the csv files to run the script on
csv_result = [i for i in glob.glob('*.{}'.format(extension))]

completed_path = r'C:\Users\U590135\.spyder-py3\projects\Soybean TLM\completed'

results_files = []
for i in csv_result:
    if 'temp_calls' in i:
        continue
    else:
        results_files.append(i)
            
for fhand in results_files:  # for each csv in the folder, run the analysis
    print('Analyzing ' + fhand + '...')

    with open(fhand, newline='') as kraken:
        reader = csv.DictReader(kraken)
        kraken_headers = reader.fieldnames
        
        panel_summaries = ['Rps1 Zygosity Call', 'Rps3 Zygosity Call', 'Rhg1 Summary Call',
                           'SCC_3 Summary Call', 'SCC_11 Summary Call']  # add any new panels summary column name to this list
    
        panel_list = []  # empty list where we can store which panels were run - for final merge sequence
        header_list = kraken_headers

        for summs in panel_summaries:
            header_list.append(summs)  # add the possible panel summaries to the headers
    
        with open('temp_calls.csv', 'w', newline='') as resultfile:
            writer = csv.DictWriter(resultfile, fieldnames=header_list)
            writer.writeheader()
    
            # try and index all of the panels to see if they are in the kraken analysis
    
            rps1_index = indexAssays('rps1_haplotypes', kraken_headers)
            rps1_len = int(len(rps1_index))
            if rps1_len > 1:
                panel_list.append('rps1_haplotypes')
    
            rps3_index = indexAssays('rps3a_haplotypes', kraken_headers)
            rps3_len = int(len(rps3_index))
            if rps3_len > 1:
                panel_list.append('rps3a_haplotypes')
    
            rhg1_index = indexAssays('rhg1_haplotypes', kraken_headers)
            rhg1_len = int(len(rhg1_index))
            if rhg1_len > 1:
                panel_list.append('rhg1_haplotypes')
    
            scc3_index = indexAssays('scc_3_haplotype', kraken_headers)
            scc3_len = int(len(scc3_index))
            if scc3_len > 1:
                panel_list.append('scc_3_haplotype')
    
            scc11_index = indexAssays('scc_11_haplotype', kraken_headers)
            scc11_len = int(len(scc11_index))
            if scc11_len > 1:
                panel_list.append('scc_11_haplotype')
                
            # example script for addition of new TLM panel below.  variables are a short form name (exe short for example)
            # and panel name must match the name exactly to the haplotype panel name at top of script
       
            exe_index = indexAssays('example: new panel name here', kraken_headers)
            exe_len = int(len(exe_index))
            if exe_len > 1:
                panel_list.append('example: new panel name here')
                
                
            # ----------------------- Reading the file and generating summary columns --------------------------------------
            
            for row in reader:
                remove_controls = ['A2', 'A02', 'A3', 'A03', 'A4', 'A04', 'A5', 'A05', 'A6', 'A06']
                if row['Well'] in remove_controls:  # remove controls from report
                    continue
                else:    
                    row['Project'] = row['Project'].encode('utf-8')  # encode the Project column as utf-8
                    

                # ------------------------  RPS1 Panel ----------------------------------
                if rps1_len > 1:
                    ''' RPS1 is present so make a dictionary out of the panel, pass into getCall to make a summary column'''
                    rps1_index.sort()  # make sure column names are in order
                    dict_rps1 = {}
    
                    for i in rps1_index:
                        dict_rps1[kraken_headers[i]] = row[kraken_headers[i]]

                    rps1_call = getCall(dict_rps1, 'rps1_haplotypes')
                    row['Rps1 Zygosity Call'] = rps1_call
    
                # --------------------------- RPS3A Panel --------------------------------
                if rps3_len > 1:
                    ''' RPS3A is present so make a dictionary out of the panel, pass into getCall to make a summary column'''
                    rps3_index.sort()  # make sure column names are in order
                    dict_rps3 = {}
    
                    for i in rps3_index:
                        dict_rps3[kraken_headers[i]] = row[kraken_headers[i]]
                    
                    rps3_call = getCall(dict_rps3, 'rps3a_haplotypes')
                    row['Rps3 Zygosity Call'] = rps3_call
    
                # ------------------------------ RHG1 Panel ---------------------------------
                if rhg1_len > 1:
                    ''' RHG1 is present so make a dictionary out of the panel, pass into getCall to make a summary column'''
                    rhg1_index.sort()  # make sure column names are in order
                    dict_rhg1 = {}
    
                    for i in rhg1_index:
                        dict_rhg1[kraken_headers[i]] = row[kraken_headers[i]]
    
                    rhg1_call = getCall(dict_rhg1, 'rhg1_haplotypes')
                    row['Rhg1 Summary Call'] = rhg1_call
    
                # ---------------------------- SCC_3 Panel --------------------------
                if scc3_len > 1:
                    ''' SCC_3 is present so convert panel into dictionary, pass into getCall to make a summary column'''
                    scc3_index.sort()  # make sure column names are in order
                    dict_scc3 = {}
    
                    for i in scc3_index:
                        dict_scc3[kraken_headers[i]] = row[kraken_headers[i]]
    
                    scc3_call = getCall(dict_scc3, 'scc_3_haplotype')
                    row['SCC_3 Summary Call'] = scc3_call
    
                # ------------------------------ SCC_11 Panel ------------------------
                if scc11_len > 1:
                    ''' SCC_11 is present so convert panel into dictionary, pass into getCall to make a summary column'''
                    scc11_index.sort()  # make sure column names are in order
                    dict_scc11 = {}
    
                    for i in scc11_index:
                        dict_scc11[kraken_headers[i]] = row[kraken_headers[i]]
    
                    scc11_call = getCall(dict_scc11, 'scc_11_haplotype')
                    row['SCC_11 Summary Call'] = scc11_call
    
                # ------------------------------- Example for Panel Addition ------------------------------
                if exe_len > 1: # this has to match your example index len variable name above
                    ''' Example is present so convert panel into dictionary, pass into getCall to make a summary column'''
                    exe_index.sort()  # this has to match index variable name above
                    dict_exe = {}  # add short-form of panel name after dict_
    
                    for i in exe_index:
                        dict_exe[kraken_headers[i]] = row[kraken_headers[i]]
                    
                 # green text in next line has to match exactly to the haplotype_panel key you add at top of script
                    exe_call = getCall(dict_exe, 'example: new panel name here') 
                    row['Example summary column name'] = exe_call
    

                writer.writerow(row)  # write new line of kraken file with summary columns now added
    
    # ----------- format the file so summary columns are after assay panels -------------------
    
    tlm_df = pd.read_csv('temp_calls.csv', encoding='UTF-8')
    # get rid of the b' prefix in project column
    tlm_df['Project'] = tlm_df['Project'].map(lambda x: x.lstrip("b'").rstrip("'"))
    
     # reformat the file and write the TLM summary report to excel
    
    krak_info = ['Box', 'Well', 'Project', 'Sample Tube Number', 'Loc Seq#', 'RowId']
    
    df_initial = tlm_df[krak_info]  # sample and well info dataframe only (no molecular data)
    
    # merge initial with rps1 dataframe if panel is in the kraken file
    
    if 'rps1_haplotypes' in panel_list:
        rps1_headers = krak_info
        for i in rps1_index:
            rps1_headers.append(kraken_headers[i])
        rps1_headers.append('Rps1 Zygosity Call')
        rps1_df = tlm_df[rps1_headers]
        merge1 = pd.merge(left=df_initial, right=rps1_df, how='left')
    else:
        tlm_df.drop(['Rps1 Zygosity Call'], axis=1, inplace=True, errors='ignore')
        merge1 = df_initial
    
    # merge that with rps3 dataframe if panel is in the kraken file
    
    if 'rps3a_haplotypes' in panel_list:
        rps3_headers = krak_info
        for i in rps3_index:
            rps3_headers.append(kraken_headers[i])
        rps3_headers.append('Rps3 Zygosity Call')
        rps3_df = tlm_df[rps3_headers]
        merge2 = pd.merge(left=merge1, right=rps3_df, how='left')
    else:
        tlm_df.drop(['Rps3 Zygosity Call'], axis=1, inplace=True, errors='ignore')
        merge2 = merge1
    
    # merge that with rhg1 dataframe if panel is in the kraken file
    
    if 'rhg1_haplotypes' in panel_list:
        rhg1_headers = krak_info
        for i in rhg1_index:
            rhg1_headers.append(kraken_headers[i])
        rhg1_headers.append('Rhg1 Summary Call')
        rhg1_df = tlm_df[rhg1_headers]
        merge3 = pd.merge(left=merge2, right=rhg1_df, how='left')
    else:
        tlm_df.drop(['Rhg1 Summary Call'], axis=1, inplace=True, errors='ignore')
        merge3 = merge2
    
    # merge that with SCC_3 dataframe if panel is in the kraken file
    
    if 'scc_3_haplotype' in panel_list:
        scc3_headers = krak_info
        for i in scc3_index:
            scc3_headers.append(kraken_headers[i])
        scc3_headers.append('SCC_3 Summary Call')
        scc3_df = tlm_df[scc3_headers]
        merge4 = pd.merge(left=merge3, right=scc3_df, how='left')
    else:
        tlm_df.drop(['SCC_3 Summary Call'], axis=1, inplace=True, errors='ignore')
        merge4 = merge3
    
    # merge that with SCC_11 spring dataframe if panel is in the kraken file
    
    if 'scc_11_haplotype' in panel_list:
        scc11_headers = krak_info
        for i in scc11_index:
            scc11_headers.append(kraken_headers[i])
        scc11_headers.append('SCC_11 Summary Call')
        scc11_df = tlm_df[scc11_headers]
        merge5 = pd.merge(left=merge4, right=scc11_df, how='left')
    else:
        tlm_df.drop(['SCC_11 Summary Call'], axis=1, inplace=True, errors='ignore')
        merge5 = merge4
        
        
    '''
    if adding another panel, you have to add it to the merge sequence here (above this paragraph).  copy the 
    section above and tweak the names of the variables to the new panel.  change Merge# by adding 1(ex. merge6).
    In merge_final below, make sure to change the 'left=' to the merge number from your new panel (ex. left=merge6)
    '''
    
    # merge with final dataframe and write the file (gets any assays not in panel tacked on to end of kraken report)
    
    merge_final = pd.merge(left=merge5, right=tlm_df, how='left')
    
    find_the_dot = fhand.find('.')
    final_name = fhand[:find_the_dot] + '_Report.xlsx'
    results_writer = ExcelWriter(final_name)
    
    os.chdir(completed_path)
    merge_final.to_excel(results_writer, 'Report', index=False)
    results_writer.save()
    
    os.chdir(path) 
    os.remove('temp_calls.csv') # delete the temporary call file
    os.remove(fhand)  # delete CSV version of kraken file
    
    print('Successfully completed {} ...\n'.format(fhand))
    
    time.sleep(2)

print('Completed all analysis')