#! /usr/bin/env python

import argparse
import sys
import pandas as pd
from datetime import date
import re 

# data type dictionaries
cov_out_data_types = {'#rname': object,
'startpos' : 'Int64',
'endpos' : 'Int64',
'numreads' : 'Int64',
'covbases' : 'Int64',
'coverage' : 'float64',
'meandepth' : 'float64',
'meanbaseseq' : 'float64',
'meanmapq' : 'float64'}


percent_cov_data_type = {'sample_name': object,
 'aligned_bases': 'Int64',
 'N_bases': 'Int64',
 'non_ambiguous_bases': 'Int64',
 'percent_coverage': 'float64'}


terra_data_table_data_types = {'index_position': 'Int64',
 'hsn': object,
 'sample_name': object,
 'sample_well': object,
 'sample_type': object,
 'index_well': object,
 'index_1': object,
 'index_2': object,
 'index_1_id': object,
 'index_2_id': object,
 'index_kit': object,
 'index_set': 'Int64',
 'plate_name': object,
 'project_name': object,
 'run_name': object,
 'library_prep': object,
 'project_type': object,
 'organism': object,
 'run_date': object,
 'read_length': 'Int64',
 'read_type': object,
 'primer_set': object,
 'platform': object,
 'instrument_id': object,
 'tag': object,
 'note': object,
 'verification_set_name': object,
 'fastq_dir': object,
 'workbook_path': object,
 'terra_data_table_path': object,
 'out_dir': object,
 'download_date': object}




#### FUNCTIONS #####
def getOptions(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument('--sample_name_array')
    parser.add_argument("--cov_out_files")
    parser.add_argument('--percent_cvg_files')
    parser.add_argument('--project_name')
    parser.add_argument('--terra_data_table_path')
   

    options = parser.parse_args(args)
    return options


def create_list_from_write_lines_input(write_lines_input):
    list = []
    with open(write_lines_input, 'r') as f:
        for line in f:
            list.append(line.strip())
    return list


def concat_cov_out(cov_out_file_list):

     # initiate dataframe for concatenation
    df = pd.DataFrame()
    sample_name_list = []
    samtools_mapped_reads_list = []
    samtools_depth_list = []
    samtools_baseq_list = []
    samtools_mapq_list = []

    #loop through bam file stats files and pull data
    for file in cov_out_file_list:
        d = pd.read_csv(file, sep = '\t', dtype = cov_out_data_types)
        if re.search('barcode', file):
            # for nanopore runs
            sample_name = re.findall('/([0-9a-zA-Z_\-\.]+)_barcode', file)[0]
        else:
            # for illumina runs
            sample_name = re.findall('/([0-9a-zA-Z_\-\.]+)_coverage.txt', file)[0]

        # pull data from samtools output
        num_reads = d.numreads[0]
        depth = d.meandepth[0]
        baseq = d.meanbaseq[0]
        mapq = d.meanmapq[0]

        sample_name_list.append(sample_name)
        samtools_mapped_reads_list.append(num_reads)
        samtools_depth_list.append(depth)
        samtools_baseq_list.append(baseq)
        samtools_mapq_list.append(mapq)

    df['sample_name'] = sample_name_list
    df['mapped_reads'] = samtools_mapped_reads_list
    df['mean_depth'] = samtools_depth_list
    df['mean_base_quality'] = samtools_baseq_list
    df['mean_map_quality'] = samtools_mapq_list

    return df

def concat_percent_cvg(percent_cvg_file_list):

    df_list = []
    for file in percent_cvg_file_list:
        d = pd.read_csv(file, dtype = percent_cov_data_type)
        df_list.append(d)

    df = pd.concat(df_list)

    return df



def concat_results(sample_name_list, terra_data_table_path, project_name, 
                   cov_out_df, percent_cvg_df):

    # set some functions for getting data formatted
    def get_sample_name_from_fasta_header(fasta_header):
        sample_name = str(re.findall('CO-CDPHE-([0-9a-zA-Z_\-\.]+)', fasta_header)[0])
        return sample_name

    def create_fasta_header(sample_name):
        return 'CO-CDPHE-%s' % sample_name
    

    # create dataframe and fill with constant strings
    df = pd.DataFrame()
    df['sample_name'] = sample_name_list
    df = df.set_index('sample_name')
    df['analysis_date'] = str(date.today())


    # read in terra_data_table
    terra_data_table = pd.read_csv(terra_data_table_path, sep = '\t', dtype = terra_data_table_data_types)
    drop_col = terra_data_table.columns.tolist()[0]
    terra_data_table = terra_data_table.drop(columns = drop_col)
    terra_data_table = terra_data_table.set_index('sample_name')


    # set index on sample_names to prepare for joining
    cov_out_df = cov_out_df.set_index('sample_name')
    percent_cvg_df = percent_cvg_df.set_index('sample_name')

    # join
    j = df.join(terra_data_table, how = 'left')
    j = j.join(percent_cvg_df, how = 'left')
    j = j.join(cov_out_df, how = 'left')
    j = j.reset_index()

    # add fasta header
    j['fasta_header'] = j.apply(lambda x:create_fasta_header(x.sample_name), axis=1)

    # add assembled column and fill in failed assembles with 0% coveage
    j.percent_coverage = j.percent_coverage.fillna(value = 0)

    def get_assembly_pass(percent_coverage):
        if percent_coverage == 0:
            return False
        if percent_coverage > 0:
            return True      
    j['assembly_pass'] = j.apply(lambda x:get_assembly_pass(x.percent_coverage), axis = 1)

    # order columns
    columns = j.columns.tolist()
    columns.sort()
    primary_columns = ['hsn', 'sample_name', 'project_name', 'plate_name', 
                       'run_name', 'analysis_date', 'run_date', 'assembly_pass', 
                       'percent_coverage']
    for column in columns:
         if column not in primary_columns:
              primary_columns.append(column)
    
    j = j[primary_columns]


    outfile = f'{project_name}_sequencing_results.csv' 
    j.to_csv(outfile, index = False)

    return j


if __name__ == '__main__':

    options = getOptions()

    sample_name_array = options.sample_name_array
    terra_data_table_path = options.terra_data_table_path
    cov_out_files = options.cov_out_files
    percent_cvg_files = options.percent_cvg_files
    project_name = options.project_name


    # create lists from the column table txt file input
    sample_name_list = create_list_from_write_lines_input(write_lines_input=sample_name_array)
    cov_out_file_list = create_list_from_write_lines_input(write_lines_input = cov_out_files)
    percent_cvg_file_list = create_list_from_write_lines_input(write_lines_input=percent_cvg_files)
    
    # concat cov_out files and percent_cvg files
    cov_out_df = concat_cov_out(cov_out_file_list=cov_out_file_list)
    percent_cvg_df = concat_percent_cvg(percent_cvg_file_list=percent_cvg_file_list)


    # create results file
    results_df = concat_results(sample_name_list = sample_name_list,
                                terra_data_table_path = terra_data_table_path,
                                project_name = project_name,
                                cov_out_df=cov_out_df,
                                percent_cvg_df=percent_cvg_df)
    
    

