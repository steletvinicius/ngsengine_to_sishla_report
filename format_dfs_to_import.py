import pandas as pd
from datetime import date

def format_df(df,client):
    cols        = ['_alelo01', '_alelo02']
    df_pivoted  = pivot_dataframe(df, cols)
    df_renamed  = rename_cols(df_pivoted, client) if (client in ['sishla', 'redome']) else df_pivoted
    df_filled   = df_replace_NAs(df_renamed)
    df_final    = append_empty_lines_to_df(df_filled) if (client in ['sishla']) else df_renamed
    
    return df_final

def pivot_dataframe(df,cols):
    # df_pivot = df.reset_index().pivot( index = 'sample_name', columns = 'locus_name', values = cols )
    df_pivot = df.reset_index().pivot( index = ['sample_name','fastq_filename','software_name','software_version'], columns = 'locus_name', values = cols )
    df_pivot.columns = ["_".join(tup) for tup in df_pivot.columns.to_flat_index()]
    return df_pivot

def rename_cols(df_pivoted, client):
    cols_df_rename_sishla = {'_alelo01_HLA-A':'01_alelo01_HLA-A_sishla',
                           '_alelo02_HLA-A':'02_alelo02_HLA-A_sishla',
                           '_alelo01_HLA-B':'03_alelo01_HLA-B_sishla',
                           '_alelo02_HLA-B':'04_alelo02_HLA-B_sishla',
                           '_alelo01_HLA-C':'05_alelo01_HLA-C_sishla',
                           '_alelo02_HLA-C':'06_alelo02_HLA-C_sishla',
                           '_alelo01_HLA-DQB1':'07_alelo01_HLA-DQB1_sishla',
                           '_alelo02_HLA-DQB1':'08_alelo02_HLA-DQB1_sishla',
						   '_alelo01_HLA-DPB1':'09_alelo01_HLA-DPB1_sishla',
                           '_alelo02_HLA-DPB1':'10_alelo02_HLA-DPB1_sishla',
                           '_alelo01_HLA-DRB1':'11_alelo01_HLA-DRB1_sishla',
                           '_alelo02_HLA-DRB1':'12_alelo02_HLA-DRB1_sishla',
                           '_alelo01_HLA-DQA1':'13_alelo01_HLA-DQA1_sishla',
                           '_alelo02_HLA-DQA1':'14_alelo02_HLA-DQA1_sishla',
                           '_alelo01_HLA-DPA1':'15_alelo01_HLA-DPA1_sishla',
                           '_alelo02_HLA-DPA1':'16_alelo02_HLA-DPA1_sishla',
                           '_alelo01_HLA-DRB3':'17_alelo01_HLA-DRB3_sishla',
                           '_alelo02_HLA-DRB3':'18_alelo02_HLA-DRB3_sishla',
                           '_alelo01_HLA-DRB4':'19_alelo01_HLA-DRB4_sishla',
                           '_alelo02_HLA-DRB4':'20_alelo02_HLA-DRB4_sishla',
                           '_alelo01_HLA-DRB5':'21_alelo01_HLA-DRB5_sishla',
                           '_alelo02_HLA-DRB5':'22_alelo02_HLA-DRB5_sishla',
                        #    'software_name': '23_software_name_sishla',
                        #     'software_version': '24_software_version_sishla',
            }
    cols_df_rename_redome = {'date':'01_date_redome',
                             '_alelo01_HLA-A':'02_alelo01_HLA-A_redome',
                           '_alelo02_HLA-A':'03_alelo02_HLA-A_redome',
                           '_alelo01_HLA-B':'04_alelo01_HLA-B_redome',
                           '_alelo02_HLA-B':'05_alelo02_HLA-B_redome',
                           '_alelo01_HLA-C':'06_alelo01_HLA-C_redome',
                           '_alelo02_HLA-C':'07_alelo02_HLA-C_redome',
                           '_alelo01_HLA-DRB1':'08_alelo01_HLA-DRB1_redome',
                           '_alelo02_HLA-DRB1':'09_alelo02_HLA-DRB1_redome',                           
                           '_alelo01_HLA-DQB1':'10_alelo01_HLA-DQB1_redome',
                           '_alelo02_HLA-DQB1':'11_alelo02_HLA-DQB1_redome',
                           '_alelo01_HLA-DRB3':'12_alelo01_HLA-DRB3_redome',
                           '_alelo02_HLA-DRB3':'13_alelo02_HLA-DRB3_redome',
                           '_alelo01_HLA-DRB4':'14_alelo01_HLA-DRB4_redome',
                           '_alelo02_HLA-DRB4':'15_alelo02_HLA-DRB4_redome',
                           '_alelo01_HLA-DRB5':'16_alelo01_HLA-DRB5_redome',
                           '_alelo02_HLA-DRB5':'17_alelo02_HLA-DRB5_redome',                           
                           '_alelo01_HLA-DQA1':'18_alelo01_HLA-DQA1_redome',
                           '_alelo02_HLA-DQA1':'19_alelo02_HLA-DQA1_redome',
                           '_alelo01_HLA-DPA1':'20_alelo01_HLA-DPA1_redome',
                           '_alelo02_HLA-DPA1':'21_alelo02_HLA-DPA1_redome',
                           '_alelo01_HLA-DPB1':'22_alelo01_HLA-DPB1_redome',
                           '_alelo02_HLA-DPB1':'23_alelo02_HLA-DPB1_redome',
            }
    #print(client)
    #cols_df_rename  = cols_df_rename_sishla
    if client == 'redome':
        cols_df_rename = cols_df_rename_redome
        #Adicionando colunas de locis ausentes
        absent_cols = [col for col in cols_df_rename.keys() if col not in df_pivoted.columns]
        df_pivoted[absent_cols] = ""
        df_pivoted["date"] = date.today()
        df_pivoted["date"]= df_pivoted["date"].apply(lambda x: x.strftime("%d/%m/%Y"))

    else:
        cols_df_rename = cols_df_rename_sishla
    # print(cols_df_rename)
    df_renamed      = df_pivoted.copy()
    #print(cols_df_rename)
    df_renamed_updt = df_renamed.rename(columns=cols_df_rename)
    #print(df_renamed_updt.columns)
    df_renamed_updt = df_renamed_updt.reindex( sorted(df_renamed_updt.columns), axis=1)
    #print(df_renamed_updt.columns)
    return df_renamed_updt


def df_replace_NAs(df):
    df.reset_index(inplace=True)
    df = df.fillna("")
    df = df.replace('nan', '')
    df = df.replace('None', '')
    return df

def append_empty_lines_to_df(df):
	df_empty = pd.DataFrame(columns = list(df.columns), index=[1,2,3,4] )
	df_empty.iloc[0:4] = ''
	df_sishla = pd.concat([df_empty, df], ignore_index = True)
	return df_sishla




