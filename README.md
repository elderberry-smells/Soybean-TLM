# Soybean Trait Linked Marker Tool

This tool will analyze the trait linked marker panels currently in the HTMA soybean molecular lab by reading a kraken export file.
The tool utilizes python 3.5.0.

version: 1.0

## Haplotype panels
there are currently 5 marker panels in the soybean molecular lab that are trait linked.  The following dictionary outlines the different 
panels and which combination of alleles for each assay will produce a 'Trait' call.  The panels themselves are made up of various markers 
that vary in length from a panel of 2 assays to a panel of 3.  

Some panels will only have 1 combination of alleles that will produce a trait, but others can have 2 combinations which are dictated in the 
list of dictionaries below.

```python
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
                    ]
    }

# for future use - look to see if panel is an unknown and return unknown
unknown_haplotypes =  {  # the following combinations will produce 'Unknown'
               'rps1_haplotypes' :[{'placeholder':'allele'}],

                'rps3a_haplotypes' : [{'placeholder':'allele'}],

                'rhg1_haplotypes' : [{'placeholder':'allele'}],

                'scc_3_haplotype' : [{'placeholder':'allele'}],
                 
                 'scc_11_haplotype' : [{'placeholder':'allele'}]
                }
                
```

### Reading a Kraken File output

Kraken exports as an xlsx file, so we have the scipt read through a specific folder and find any xlsx files that need analysis.  
It converts the files to a csv using pandas, and then makes a list of csv files for running through the TLM tool analysis for 
a report output.

```python
import csv
import pandas as pd
from pandas import ExcelWriter
import glob
import os
import time

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
```

### Running the files through the TLM tool

The file is read for the fieldnames as a header list and the possible panel summary headers are appended to that list for a 
temporary file creation later.  Each file that is in the folder for analysis will go through the following steps in sequence.

```python
for fhand in results_files:  # for each csv in the folder, run the analysis
    print('Analyzing ' + fhand + '...')

    with open(fhand, newline='') as kraken:
        reader = csv.DictReader(kraken)
        kraken_headers = reader.fieldnames
        
        panel_summaries = ['Rps1 Zygosity Call', 'Rps3 Zygosity Call', 'Rhg1 Summary',
                           'SCC_3 Summary', 'SCC_11 Summary']  
    
        panel_list = []  # empty list where we can store which panels were run - for final merge sequence
        header_list = kraken_headers

        for summs in panel_summaries:
            header_list.append(summs)  # add the possible panel summaries to the headers
        
        with open('temp_calls.csv', 'w', newline='') as resultfile:
            writer = csv.DictWriter(resultfile, fieldnames=header_list)
            writer.writeheader()
```

### Index the kraken headers for the assay panels

Run the kraken_header list through the indexAssays function to find the location of **every TLM panel present** in the kraken file.  
The function will return the column location (index) where the panel is present for all assays.

```python
# example: try and index rps1_haplotype panel from header lists and return a list of column locations
kraken_headers = ['Box', 'Well', 'Project', 'Sample Tube Number', 'Loc Seq#',
                  'DAS_Gm03_4458273_G_A Zygosity Call', 'DAS_Gm03_3838302_A_T Zygosity Call']
panel_list =[]

rps1_index = indexAssays('rps1_haplotypes', kraken_headers)
rps1_len = int(len(rps1_index))
if rps1_len > 1:
    panel_list.append('rps1_haplotypes')
```

```python
def indexAssays(panel, headers):
    
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

print(rsp1_index)
```    

### removing the control wells and converting the 'Project' well to utf-8

```python
for row in reader:
                remove_controls = ['A2', 'A02', 'A3', 'A03', 'A4', 'A04', 'A5', 'A05', 'A6', 'A06']
                if row['Well'] in remove_controls:  # remove controls from report
                    continue
                else:    
                    row['Project'] = row['Project'].encode('utf-8')  # encode the Project column as utf-8
```

### The getCall function
this function dictates the trait linked call of the haplotypes in the panels present in the panel.  There are a few call varieties 
that can be a result:  Trait, Wildtype, Seg, No call, No data, and blank.

The data that gets passed into the getCall function is a dictionary that is created by isolating the panels into a dictionary (from 
the previous indexing step).

```python
def getCall(dic, panel):
    '''determine the samples call for each panel in kraken file'''
    global haplotype_panel
    global unknown_haplotypes
     
    seg = 0
    dupe = 0
    no_data = 0
    empty = 0
    
    #determine if the combination of assays is the Trait haplotype
    for hap in haplotype_panel[panel]:
        if dic == hap:
            return 'Trait'
                 
    # for future use: determine if the combination of assays is the unknown haplotype
    for hap in unknown_haplotypes[panel]:
        if dic == hap:
            return 'Unknown'                         
                
    # determine if there is a seg/no call or no data in the panel
    for assay, allele in dic.items():
        if allele == 'No Data':
            no_data += 1
        elif allele == 'Empty':
            empty += 1
        elif allele == 'Unused':
            no_data += 1
        elif allele == 'Dupe':
            dupe += 1
        elif allele[0] != allele[-1]: 
            seg += 1
             
    # using the counts of each type, return the call        
    if empty > 0:
        return ''
    elif no_data > 0:
        return 'No Data'
    elif dupe > 0:
        return 'No Call'
    elif seg > 0:
        return 'Seg'
    else:
        return 'Wildtype'
        
# Example of different haplotypes passing into the getCall function

example1 = {'DAS_Gm03_4458273_G_A Zygosity Call': 'G:G', 'DAS_Gm03_3838302_A_T Zygosity Call': 'A:A'}
example2 = {'DAS_Gm03_4458273_G_A Zygosity Call': 'A:A', 'DAS_Gm03_3838302_A_T Zygosity Call': 'T:T'}
example3 = {'DAS_Gm03_4458273_G_A Zygosity Call': 'A:A', 'DAS_Gm03_3838302_A_T Zygosity Call': 'A:A'}
example4 = {'DAS_Gm03_4458273_G_A Zygosity Call': 'G:A', 'DAS_Gm03_3838302_A_T Zygosity Call': 'A:A'}
example5 = {'DAS_Gm03_4458273_G_A Zygosity Call': 'G:G', 'DAS_Gm03_3838302_A_T Zygosity Call': 'Dupe'}

example1_call = getCall(example1, 'rps1_haplotypes')
example2_call = getCall(example2, 'rps1_haplotypes')
example3_call = getCall(example3, 'rps1_haplotypes')
example4_call = getCall(example4, 'rps1_haplotypes')
example5_call = getCall(example5, 'rps1_haplotypes')

print(example1_call, example2_call, example3_call, example4_call, example5_call)
```

The various panels will pass similar to the example above but use the data from the csv being read to populate the dictionary 
being passed into the getCall function.  They are all referencing the specific panel key in the haplotype_panel to determine 
if the Trait haplotype is present.  If not, it determines if there is any no data/seg/dupe and anything left is a wildtype.

```python
if rps1_len > 1:
    ''' RPS1 is present so make a dictionary out of the panel, pass into getCall to make a summary column'''
    rps1_index.sort()  # make sure column names are in order
    dict_rps1 = {}

    for i in rps1_index:
        dict_rps1[kraken_headers[i]] = row[kraken_headers[i]]

    rps1_call = getCall(dict_rps1, 'rps1_haplotypes')
    row['Rps1 Zygosity Call'] = rps1_call

if rps3_len > 1:
    ''' RPS3A is present so make a dictionary out of the panel, pass into getCall to make a summary column'''
    rps3_index.sort()  # make sure column names are in order
    dict_rps3 = {}

    for i in rps3_index:
        dict_rps3[kraken_headers[i]] = row[kraken_headers[i]]

    rps3_call = getCall(dict_rps3, 'rps3a_haplotypes')
    row['Rps3 Zygosity Call'] = rps3_call
```

## Merging the data together
The final report will be made using pandas dataframes of each individual haplotype present.  This is to ensure the summary column 
follows the panel of assays rather than at the end of the report.

```python
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

# repeat process for each panel in sequence until you have gone through all panels (adding 1 to merge number each time). 
```

## final merge sequence and writing of new report 
the merging is completed by merging the data that was not a part of the TLM panels (regular zygosity data).  Once that is in a dataframe, 
write the file as xlsx to a "completed" folder, and remove any copies of the file from the initial folder.

```python
# merge with final dataframe and write the file (gets any assays not in panel tacked on to end of kraken report)
    
    merge_final = pd.merge(left=merge5, right=tlm_df, how='left')
    
    find_the_dot = fhand.find('.')
    final_name = fhand[:find_the_dot] + '.xlsx'
    results_writer = ExcelWriter(final_name)
    
    os.chdir(completed_path)
    merge_final.to_excel(results_writer, 'Report', index=False)
    results_writer.save()
    
    os.chdir(path) 
    os.remove('temp_calls.csv') # delete the temporary call file
    os.remove(fhand)
```
