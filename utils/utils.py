import os
import pandas as pd
from pathlib import Path

def check_dir(dir):
    if os.path.exists(dir):
        pass
    else:
        os.makedirs(dir,exist_ok=True)
    

def check_file(file,delete=False):
    if os.path.exists(file) and delete==True:
        os.remove(file)
        return False
    elif os.path.exists(file) and delete==False:
        return True
    else:
        return False


def handle_pos_name(pos_name):
    s_item=pos_name.strip().split("_")
    chrom=str(s_item[0])
    pos=int(s_item[1])
    ref=str(s_item[2])
    alt=str(s_item[3])

    return chrom,pos,ref,alt


def str2dict(q):
    """Convert quality string to dictionary"""
    dict = {int(i.split(':')[0]):int(i.split(':')[1]) for i in q.split(',') if q !="NA"}
    return dict
    

def str2bool(v):
    """ensure boolean input"""
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ValueError('Boolean value expected.')

