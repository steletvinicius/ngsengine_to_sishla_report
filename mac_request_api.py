import requests

# # API Documentation webpage: https://hml.nmdp.org/mac/raml/allele-code.raml
# # https://github.com/nmdp-bioinformatics/multiple-allele-code/tree/master/client/bash

def encode_mac(string, imgt_version):
    print(f'encoding {string}')
    api_url  = "https://hml.nmdp.org/mac/api/"
    email = "labhlainca@gmail.com"
    imgtversion = imgt_version
    trial_status = "false"

    encode_url = f'{api_url}encode?trialRun={trial_status}&email={email}&imgtHlaRelease={imgtversion}'
    url = api_url
    response = requests.post(encode_url, data=string)
    typing = response.text
    return typing
    #return f'{typing}_encoded'
# Testando se MAC consegue lidar com GLStrings
# imgt_version = '3.54.0'
# print(encode_mac("DPB1*85:01:01:01/DPB1*85:01:01:02+DPB1*105:01:01:01/DPB1*105:01:01:02/DPB1*105:01:01:03/DPB1*105:01:01:04/DPB1*105:01:01:05/DPB1*105:01:01:06/DPB1*105:01:01:07/DPB1*105:01:01:08/DPB1*105:01:01:09/DPB1*105:01:01:10/DPB1*105:01:01:11/DPB1*105:01:01:12/DPB1*105:01:01:13/DPB1*105:01:01:14/DPB1*105:01:01:15/DPB1*105:01:01:16/DPB1*105:01:01:17/DPB1*105:01:01:18/DPB1*665:01:01/DPB1*1072:01|DPB1*463:01:01:01/DPB1*463:01:01:02/DPB1*463:01:01:03+DPB1*901:01",imgt_version))

