# https://github.com/nmdp-bioinformatics/py-ard#install-from-pypi
# pip install py-ard

import pyard
import re
from mac_request_api import *

# Initialize ARD object with a version of IMGT HLA database
def init_pyard(imgt_version):
  return pyard.init(imgt_version.replace(".",""), data_dir='./py_ard_db')

def null_alleles_to_investigate():
    null_list = [
        "A*01:04N","A*01:123N","A*01:15N","A*01:16N","A*01:57N","A*02:113N","A*02:125N",
        "A*02:227N","A*02:514N","A*02:53N","A*02:83N","A*02:94N","A*03:21N","A*11:109N",
        "A*11:21N","A*23:11N","A*23:19N","A*24:09N","A*24:11N","A*24:252N","A*24:36N",
        "A*24:84N","A*24:90N","A*25:12N","A*30:70N","A*30:78N","A*31:14N","A*31:60N",
        "A*32:27N","A*32:45N","A*34:10N","A*68:18N",
        "B*07:181N","B*07:67N","B*14:07N","B*15:01:01:02N","B*15:181N","B*15:190N",
        "B*15:79N","B*35:165N","B*37:03N","B*37:42N","B*39:40N","B*40:142N","B*40:155N",
        "B*40:22N","B*44:23N","B*51:11N",
        "C*02:38N","C*02:92N","C*04:09N","C*04:93N","C*04:95N","C*05:07N","C*05:99N",
        "C*06:16N","C*06:79N","C*07:104N","C*07:198N","C*07:227N","C*07:32N","C*07:33N",
        "C*07:452N","C*07:55N","C*07:61N","C*08:127N","C*15:122N","C*16:30N",
        "DPB1*120:01N","DPB1*154:01N","DPB1*161:01N","DPB1*218:01N","DPB1*357:01N",
        "DPB1*570:01N","DPB1*61:01N","DPB1*64:01N","DQB1*02:18N","DQB1*02:20N","DQB1*03:118N",
        "DQB1*06:144N","DQB1*06:26N","DQB1*06:75N","DQB1*06:77N","DRB1*07:10N","DRB1*07:26N",
        "DRB1*12:24N","DRB4*01:03:01:02N","DRB4*01:16N","DRB4*02:01N","DRB4*03:01N","DRB5*01:08N","DRB5*01:10N",
    ]
    return null_list

def reduce_ard(ard, gl_string, option):
  #print(gl_string)
  return ard.redux(gl_string, option)

def typing_3_fields(ard, gl_string):
    return reduce_ard(ard, gl_string,'exon')

def typing_2_fields(ard, gl_string):
    return reduce_ard(ard, gl_string,'U2')

def typing_p_group(ard, gl_string):
    return reduce_ard(ard, gl_string,'P')

""" 
def filter_typing(ard, gl_string):
    null_list = null_alleles_to_investigate()
    #discard None results
    if gl_string == None:
        return ""
    else:
        if ((gl_string.find("/") == -1) and (gl_string.find("|") == -1) and (gl_string.find("|") == "N")):
            return typing_3_fields(ard, gl_string)
        else:
            # Dealing with phasing ambiguities (changed in 22/05/2024 due to )
            if gl_string.find("|") != -1:
                strings1 = []
                strings2 = []
                for st in gl_string.split("|"):
                    strings1.append(st.split('+')[0]),strings2.append(st.split('+')[1])
                alleles1=[]
                alleles2=[]

                for string in strings1:
                    alleles = (string.split('/'))
                    alleles1 += alleles

                for string in strings2:
                    alleles = (string.split('/'))
                    alleles2 += alleles
                    filtered_typing1 = [allele for allele in alleles1 if (allele[-1] != "N") or (allele in null_list)]
                    filtered_typing2 = [allele for allele in alleles2 if (allele[-1] != "N") or (allele in null_list)]    
                    typing = f'{"/".join(filtered_typing1)}+{"/".join(filtered_typing2)}'
            
            # Dealing with other ambiguities
            else:
                if gl_string.find("+") != -1:
                    string1, string2 = gl_string.split("+")
                    alleles1, alleles2   = string1.split('/'), string2.split('/')
                    filtered_typing1 = [allele for allele in alleles1 if (allele[-1] != "N") or (allele in null_list)]
                    filtered_typing2 = [allele for allele in alleles2 if (allele[-1] != "N") or (allele in null_list)]    
                    typing = f'{"/".join(filtered_typing1)}+{"/".join(filtered_typing2)}'
            
                #DRB345
                else:
                    alleles1 = gl_string.split('/')
                    filtered_typing1 = [allele for allele in alleles1 if (allele[-1] != "N") or (allele in null_list)]
                    typing = "/".join(filtered_typing1)
            
            typing = typing_3_fields(ard, typing)
            if (typing.find("/") != -1):
                if typing.find("+"):
                    allele1,allele2 = typing.split("+")
                    allele1 = typing_2_fields(ard, allele1) if allele1.find("/") != -1 else allele1
                    allele2 = typing_2_fields(ard, allele2) if allele2.find("/") != -1 else allele2
                    typing = f'{allele1}+{allele2}'

                else:
                    typing = typing_2_fields(ard, typing)

            return typing
 """        

""" 
def cemo_typing(ard, gl_string):
    typing = filter_typing(ard, gl_string)
    if (typing.find("/") != -1):
        if typing.find("+"):
            allele1,allele2 = typing.split("+")
            allele1 = typing_p_group(ard, allele1) if allele1.find("/") != -1 else allele1
            allele2 = typing_p_group(ard, allele2) if allele2.find("/") != -1 else allele2
            typing_p = f'{allele1}+{allele2}'

        else:
            typing_p = typing_p_group(ard, typing)
        
            if typing_p.find("/") !=-1:
            return typing_p + "Resolver-ambiguidade-grupo-P"
        else:
        return typing_p
        
    return typing """

# Versao antiga substituida pela funcao que usa exclusivamente o mac request
""" 
def redome_typing(ard, gl_string, imgt_version):
    typing = filter_typing(ard, gl_string)
    if (typing.find("/") != -1):
        if typing.find("+"):
            allele1,allele2 = typing.split("+")
            allele1 = typing_p_group(ard, allele1) if allele1.find("/") != -1 else allele1
            allele2 = typing_p_group(ard, allele2) if allele2.find("/") != -1 else allele2
            typing_p = f'{allele1}+{allele2}'
        else:
            typing_p = typing_p_group(ard, typing)
        
            if typing_p.find("/") !=-1:
            return typing_p + "Resolver-ambiguidade-grupo-P"
        else:
        return encode_mac(typing, imgt_version)
    return typing
 """

""" 
def hatk_format(ard, gl_string):
    if gl_string == None:
        return None
    else:
        if gl_string.find("|") != -1:
            gl_string = gl_string.split("|")[0]
        typing = reduce_ard(ard, gl_string, 'U2')
        if typing.find("/") != -1:
            typing_01 = typing.split("+")[0].split("/")[0]
            typing_02 = typing.split("+")[1].split("/")[0]
            typing = f'{typing_01}+{typing_02}'
    return typing        
 """
#imgt_version = '3.54.0'
#print(imgt_version)
#ard = init_pyard(imgt_version)
#gl_string = "DPB1*85:01:01:01/DPB1*85:01:01:02+DPB1*105:01:01:01/DPB1*105:01:01:02/DPB1*105:01:01:03/DPB1*105:01:01:04/DPB1*105:01:01:05/DPB1*105:01:01:06/DPB1*105:01:01:07/DPB1*105:01:01:08/DPB1*105:01:01:09/DPB1*105:01:01:10/DPB1*105:01:01:11/DPB1*105:01:01:12/DPB1*105:01:01:13/DPB1*105:01:01:14/DPB1*105:01:01:15/DPB1*105:01:01:16/DPB1*105:01:01:17/DPB1*105:01:01:18/DPB1*665:01:01/DPB1*1072:01|DPB1*463:01:01:01/DPB1*463:01:01:02/DPB1*463:01:01:03+DPB1*901:01"

#print(typing_3_fields(ard, gl_string))
#print(typing_2_fields(ard, gl_string))
