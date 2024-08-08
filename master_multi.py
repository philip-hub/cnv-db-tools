"""
Make a dataframe with the following

arm, cn, ai, sample

"""
#imports
import pandas as pd; import numpy as np; import re
import os


folder = "karfiles/"

sample_number = 0

def getAICN(kar_data):
    ai = kar_data['ai']
    dcn = kar_data['cn']  
    arms_kar_data = kar_data['arm']
    
    ai_cn_data = pd.DataFrame({
        'arm':arms_kar_data,
        'cn':dcn,
        'ai':ai
    })
    
    
    return ai_cn_data



