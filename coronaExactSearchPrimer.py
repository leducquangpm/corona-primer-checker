from Bio.Seq import Seq
from Bio import SeqIO


import time

import pandas as pd
from numpy import *
import re
import itertools
import csv
import sys
import os, shutil, glob
import argparse
def count_gap(text):
    if text.strip() == "": # To take of care of all space input
        return 0
    #count = sum([1 if char=='N' else 0 for char in text ])
    count = text.count('N') # Your error is here, Only check for 1 space instead of 3 spaces
    count = count+text.count('n')
    total_chars = len(text)

    return round(count / total_chars,3)*100
def ReverseComplement(Pattern):
    str=""
    for c in Pattern:
        if c=="A":
            str=str+"T"
        elif c=="T":
            str=str+"A"
        elif c=="G":
            str=str+"C"
        elif c=="C":
            str=str+"G"
        else:
            str=str+c
    return str[::-1]

def check_pos(list_primer_pos,ref_d):
    #check order
    for i in range(1,len(list_primer_pos)):
        if list_primer_pos[i]-list_primer_pos[i-1]<0:
            return False,'wrong order'
    #check length
    d=list_primer_pos[len(list_primer_pos)-1]-list_primer_pos[0]
    if d<ref_d-50 or d>ref_d+50:
        return False,'wrong distance'
    return True,'ok'

def pipeline(args):
    #read primer text:
    primers=[]
    file1 = open(args.primer, 'r')
    num_primer = int(file1.readline())
    # Strips the newline character
    for i in range(num_primer):
        primer={}
        primer['name']=file1.readline().strip()
        #print(primer['name'])
        tok=file1.readline().strip().split(' ')
        #print(tok)
        length =int(tok[0])
        num_p=int(tok[1])
        primer['length']=length
        primer['primer']=[]
        for j in range(num_p):
            p={}
            tok=file1.readline().strip().split(' ')
            p['order']=int(tok[0])
            p['type']=tok[1].strip()
            p['seq']=tok[2].strip().upper()
            if p['order']<0:
                p['seq']=ReverseComplement(p['seq'])
            primer['primer'].append(p)
        primers.append(primer)



    #ret=parseHandler(r'G:\Downloads\primer_blast.xml')
    #print(ret)
    #simple check
    list_sample=[]
    most_gap=0
    for record in SeqIO.parse(args.db, "fasta"):
        sample={}
        sample['title']=record.description
        sample['seq']=record.seq.upper()
        #filter gap
        perc_gap=count_gap(sample['seq'])
        #print(perc_gap)
        if perc_gap>most_gap:
            most_gap=perc_gap
        if perc_gap<=1 and len(record.seq)>20000:
            list_sample.append(sample)
    print(len(list_sample))
    print(most_gap)
    count_mm={}
    set_mm={}
    file=open(os.path.join(args.output,'report_FN.txt'),'w')
    for group in primers:
        count_match=0
        count_not_match=0
        count_not_contain=0

        set_mm[group['name']]=set()
        for sample in list_sample:
            list_post=[]
            isContain=True
            isMatch=False
            if not sample['title'] in count_mm:
                count_mm[sample['title']]={}
                count_mm[sample['title']]['count']=0
                count_mm[sample['title']]['primer']=''
            primer_miss=set()
            for i in range(len(group['primer'])):
                list_p=[m.start() for m in re.finditer('(?='+str(group['primer'][i]['seq'])+')', str(sample['seq']))]

                if len(list_p)<=0:
                    primer_miss.add(group['primer'][i]['order'])
                    isContain=False
                list_post.append(list_p)

            if isContain:
                for p in itertools.product(*list_post):
                    ret_c,mes=check_pos(p,group['length'])
                    if ret_c :
                        count_match=count_match+1
                        isMatch=True
                        break
            if not isContain:
                count_not_contain=count_not_contain+1
                #print(group['name']+" not contain in "+sample['title'])
                file.write(group['name']+'\t miss:'+str(primer_miss)+"\t"+sample['title']+"\n")
                count_mm[sample['title']]['count']=count_mm[sample['title']]['count']+1
                count_mm[sample['title']]['primer']=count_mm[sample['title']]['primer']+','+group['name']
                set_mm[group['name']].add(sample['title'])
            if isContain and not isMatch:
                count_not_match=count_not_match+1
                #print(group['name']+" not match in position in "+sample['title'])
                file.write(group['name']+'\t wrong order \t'+sample['title']+"\n")
                count_mm[sample['title']]['count']=count_mm[sample['title']]['count']+1
                count_mm[sample['title']]['primer']=count_mm[sample['title']]['primer']+','+group['name']
                set_mm[group['name']].add(sample['title'])
        print('group '+group['name']+':'+str(count_match)+","+str(count_not_contain)+","+str(count_not_match))
    file.close()

    max_v=0
    title_max=''
    primer_miss=''
    for key in count_mm:
        if count_mm[key]['count']>max_v:
            max_v=count_mm[key]['count']
            title_max=key
            primer_miss=count_mm[key]['primer']
    print(title_max+":"+str(max_v))
    print(primer_miss)
    set_E=set()
    set_N1=set()
    set_N2=set()
    set_N3=set()
    set_N4=set()
    set_Rd1=set()
    set_Rd2=set()
    for key in set_mm:
        if 'E' in key:
            if len(set_E)==0:
                set_E=set_mm[key]
            else:
                set_E=set_E&set_mm[key]
        if 'N1' in key:
            if len(set_N1)==0:
                set_N1=set_mm[key]
            else:
                set_N1=set_N1&set_mm[key]
        if 'N2' in key:
            if len(set_N2)==0:
                set_N2=set_mm[key]
            else:
                set_N2=set_N2&set_mm[key]
        if 'N3' in key:
            if len(set_N3)==0:
                set_N3=set_mm[key]
            else:
                set_N3=set_N3&set_mm[key]
        if 'N4' in key:
            if len(set_N4)==0:
                set_N4=set_mm[key]
            else:
                set_N4=set_N4&set_mm[key]
        if 'Rd1' in key:
            if len(set_Rd1)==0:
                set_Rd1=set_mm[key]
            else:
                set_Rd1=set_Rd1&set_mm[key]
        if 'Rd2' in key:
            if len(set_Rd2)==0:
                set_Rd2=set_mm[key]
            else:
                set_Rd2=set_Rd2&set_mm[key]
    bag_miss={'E':set_E,'N1':set_N1,'N2':set_N2,'N3':set_N3,'N4':set_N4,'Rd1':set_Rd1,'Rd2':set_Rd2}
    f=open(os.path.join(args.output,'report_FN_exact_by_group.tab'),'w')
    for k in bag_miss:
        for s in bag_miss[k]:
            f.write(k+'\t'+s+'\n')
    f.close
def main(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        prog='corona_fn_checker',
        description='Tool for checking False Negative error when testing specified primers with corona')
    subparsers = parser.add_subparsers(title='sub command', help='sub command help')
    run_cmd = subparsers.add_parser(
        'run', description='Check haplotype status with specified primer', help='Check FN',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    run_cmd.set_defaults(func=pipeline)
    run_cmd.add_argument('--primer', help='Primer list in text file',type=str)
    run_cmd.add_argument('--db', help='Corona db file',type=str)
    run_cmd.add_argument('--output', help='Output folder', type=str)
    args = parser.parse_args(arguments)
    return args.func(args)
if __name__ == "__main__":
    main()
