from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW
import time
from Bio.Blast import NCBIXML

from numpy import *
import re
import itertools
import csv
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
    
def parseHandler(file_xml):
    dict_ret={}
    result_handle = open(file_xml)
    blast_records = NCBIXML.parse(result_handle)
    blast_records=list(blast_records)
    #print(len(blast_records))
   
    for blast_record in blast_records:
        list_ret={}
        #for desc in blast_record.descriptions:
        #    print(desc.title)
        tok=blast_record.query.split('|')
        name=tok[0].strip()
        prime_index=tok[1].strip()
        if not blast_record.query in dict_ret:
            dict_ret[name]={}
        if not prime_index in dict_ret[name]:
            dict_ret[name][prime_index]={}
        for align in blast_record.alignments:
            #print(align.title) 
            dic={}
            dic['title']=align.title
            m=align.title.split('|')
            for i in range(len(m)):
                if m[i]=='emb':
                    dic['id']=m[i+1]
                if m[i]=='ref':
                    dic['id']=m[i+1]
                if m[i]=='dbj':
                    dic['id']=m[i+1]
                if m[i]=='gi':
                    dic['id']=m[i+1]
                if m[i]=='gb':
                    dic['id']=m[i+1]
            dic['hsps']=[]
            overThreashold=False
            for hsp in align.hsps:
                d={}
                d['score']=hsp.score
                d['expect']=hsp.expect
                d['indentity']=hsp.identities/blast_record.query_length
                if d['indentity']>=1:
                    overThreashold=True
                d['subject_start']=hsp.sbjct_start
                d['subject_end']=hsp.sbjct_start+blast_record.query_length
                dic['hsps'].append(d)
            if overThreashold:
                list_ret[dic['id']]=dic
        dict_ret[name][prime_index]['sample']=list_ret
    return dict_ret

#read primer text:
primers=[]
file1 = open('primerUltramp.txt', 'r') 
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
f=open("primer.fasta",'w')
for p in primers:
    for primer in p['primer']:
        f.write('>'+p['name']+"|"+str(primer['order'])+"\n")
        f.write(primer['seq']+"\n")
f.close()

#ret=parseHandler(r'G:\Downloads\primer_blast.xml')
#print(ret)
#simple check
list_sample=[]
most_gap=0
for record in SeqIO.parse("SarsCov2_NCBI_785ncbi.fasta", "fasta"):
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
file=open('report_FN.txt','w')
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
f=open('report_FN_exact_by_group.tab','w')
for k in bag_miss:
    for s in bag_miss[k]:
        f.write(k+'\t'+s+'\n')
f.close 
