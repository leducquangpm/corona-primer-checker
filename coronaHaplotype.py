import subprocess
import os, shutil, glob
import re
import sys
import itertools
import csv
import argparse
import json
from Bio import SeqIO
def count_gap(text):
    if text.strip() == "": # To take of care of all space input
        return 0
    #count = sum([1 if char=='N' else 0 for char in text ])
    count = text.count('N') # Your error is here, Only check for 1 space instead of 3 spaces
    count =count+ text.count('n')
    total_chars = len(text)

    return round(count / total_chars,3)*100
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
def getGoodCombination(pos_start,pos_end,length):
    list_pos=[]
    list_pos.append(pos_start)
    list_pos.append(pos_end)
    correct_p=[]
    for p in itertools.product(*list_pos):
        retc,mes=check_pos(p,length)
        if retc:
            correct_p=p
            break
    return correct_p
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
        elif c=='S':
            str=str+"S"
        elif c=='R':
            str=str+"Y"
        elif c=='Y':
            str=str+"R"
        elif c=='K':
            str=str+"M"
        elif c=='M':
            str=str+"K"
        else:
            str=str+c
    return str[::-1]
def combileFastaFile(folder_in,file_out):
    list_seq=[]

    for root, dirs, files in os.walk(folder_in):
        for _file in files:
            if _file.endswith(('.fasta')):
                #print(str(root)+'/'+_file)
                for seq in SeqIO.parse(str(root)+'/'+_file,'fasta'):
                    perc_gap=count_gap(seq.seq)
                    if perc_gap<=0 and len(seq.seq)>20000:
                        list_seq.append(seq)
        SeqIO.write(list_seq,file_out,'fasta')
def filterFastaFile(fasta_in,file_out):
    list_seq=[]
    for seq in SeqIO.parse(fasta_in,'fasta'):
        perc_gap=count_gap(seq.seq)
        if perc_gap<=1 and len(seq.seq)>20000:
            
            list_seq.append(seq)
    SeqIO.write(list_seq,file_out,'fasta')
    return file_out
def getAMP(region,primerFi,primerRo):
    s_pos,mm,m=FittingAlignment(region,primerFi,1,1)
    e_pos,mm,m=FittingAlignment(region,primerRo,1,1)
    e_pos=e_pos+len(primerRo)
    return region[s_pos:e_pos],s_pos,e_pos
def setupdb(db_file):
    """
    make blast database from fasta file in db folder,
    :param : fasta file (with folder'folder is the name of db and filename is 'sequences')
    :return:
    """
    #seqfile='/mnt/data/coronacheck/sarscov2.fasta'
    #gisaid_dir='/home/quang/Downloads/gisaid4000'
    #db_file='/mnt/data/coronacheck/SarsCov2_NCBI_785ncbi.fasta'
    #FN_ref_dir='/mnt/data/coronacheck/fpsamples'

    list_seq=[]

    # for root, dirs, files in os.walk(FN_ref_dir):
    #     for _file in files:
    #         if _file.endswith(('.fasta')):
    #             #print(str(root)+'/'+_file)
    #             for seq in SeqIO.parse(str(root)+'/'+_file,'fasta'):
    #                 perc_gap=count_gap(seq.seq)
    #                 if perc_gap<=0 and len(seq.seq)>20000:
    #                    list_seq.append(seq)
    # SeqIO.write(list_seq,db_file,'fasta')
    # print(len(list_seq))
    for seq in SeqIO.parse(db_file,'fasta'):
        perc_gap=count_gap(seq.seq)
        if perc_gap<=1 and len(seq.seq)>20000:
            list_seq.append(seq)
    print(len(list_seq))
    SeqIO.write(list_seq,db_file,'fasta')
    cmd="makeblastdb -in {path} -title {name} -dbtype {type} -logfile /dev/null".format(
                            path=db_file,
                            name='corona',
                            type='nucl'

    )
    print (cmd)
    os.system(cmd)
    return db_file
def setupdbRef(ref_file):
    cmd="makeblastdb -in {path} -title {name} -dbtype {type} -logfile /dev/null".format(
                            path=ref_file,
                            name='corona',
                            type='nucl'

    )
    print (cmd)
    os.system(cmd)
import uuid 
def blast(sample,db, identity=70, threads=10, mincov=70,dbtype='nucl'):
    """
    Call blastn with params
    :param query_file (in fasta), db (blast indexed db), number of threads and identity
    :return: list BLASTFields objects
    """
    #check db is indexed
    #dbfile=os.path.join(db_folder, 'sequences')
    # run blastn
    uid=str(uuid.uuid1().hex)
    temp_fasta=uid+'_temp.fasta'
    temp_tab=uid+'_temp.tab'
    f=open(temp_fasta,'w')
    f.write('>primer\n')
    f.write(sample)
    f.close()
    cmd ='blastn -query {query} -task blastn-short -max_target_seqs 1000000 -perc_identity {identity} -db {db} -outfmt \'6 qseqid qstart qend qlen sseqid sstart send slen sstrand length pident mismatch qseq sseq stitle bitscore\'  -num_threads {threads} >{output}'.format(

        query=temp_fasta,
        identity=identity,
        db=db,
        threads=threads,
        output=temp_tab

    )
    print(cmd)
    os.system(cmd)
    #parse result
    f=open(temp_tab)
    line = f.readline()
    result=[]
    while line:
        #result.append(line)
        t=line.strip().split('\t')
        blast_fields={'qseqid':t[0], 'qstart':t[1], 'qend':t[2], 'qlen':t[3],\
         'sseqid':t[4], 'sstart':t[5], 'send':t[6], 'slen':t[7], 'sstrand':t[8],\
           'length':t[9],'pident':t[10], 'mismatch':t[11], 'qseq':t[12], 'sseq':t[13],\
           'stitle':t[14],'bitscore':t[15]}
        result.append(blast_fields)
        line = f.readline()
    f.close()
    if os.path.exists(temp_tab):
        os.remove(temp_tab)
    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)
    return result
def scoreTable(a,b):
    if a==b:
        return 1
    else:
        return -1
def FittingAlignment(v,w,score,sigma):
    if v=='' or w=='':
        return 0,'',0
    s=[[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    bk=[[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    s[0][0]=0
    for i in range(1,len(v)+1):
        s[i][0]=0
    for j in range(1,len(w)+1):
        s[0][j]=s[0][j-1]-1
    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
            if v[i-1]==w[j-1]:
                s[i][j]=max(s[i-1][j]-sigma,s[i][j-1]-sigma,s[i-1][j-1]+score)
            else:
                s[i][j]=max(s[i-1][j]-sigma,s[i][j-1]-sigma,s[i-1][j-1]-sigma)
            if s[i][j]==s[i-1][j]-sigma:
                bk[i][j]=1
            elif s[i][j]==s[i][j-1]-sigma:
                bk[i][j]=2
            else:
                bk[i][j]=3
    AlignV=''
    AlignW=''
    mi=0
    mj=len(w)
    for i in range(0,len(v)+1):
        if s[i][len(w)] >s[mi][len(w)]:
            mi=i
    i=mi
    j=mj
    while True:
        if j>0:
            if bk[i][j] ==3:

                AlignV=v[i-1]+AlignV
                AlignW=w[j-1]+AlignW
                i=i-1
                j=j-1
            elif bk[i][j] ==1:
            
                AlignW='-'+AlignW
                AlignV=v[i-1]+AlignV
                i =i-1
            else:
            
                AlignW=w[j-1]+AlignW
                AlignV='-'+AlignV
                j =j-1
        
        else:
            break
          
    mi=i        
    mm=''
    m=0
    for i in range(len(AlignW)):
        if not AlignV[i]==AlignW[i]:
            mm=mm+str(i+1)+'('+AlignV[i]+'->'+AlignW[i]+')'+';'
        else:
            m=m+1
    
    count_gap=0
    for i in range(len(AlignV)):
        if not AlignV[i]=='-':
            break
        else:
            count_gap=count_gap+1
    mi=mi-count_gap
    return mi,mm,m/(len(w))
def GlobalAlignment(v,w,score,sigma):
    if v=='' or w=='':
        return '',0
    s=[[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    bk=[[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
    s[0][0]=0
    for i in range(1,len(v)+1):
        s[i][0]=s[i-1][0]-sigma
    for j in range(1,len(w)+1):
        s[0][j]=s[0][j-1]-sigma
    for i in range(1,len(v)+1):
        for j in range(1,len(w)+1):
            # s[i][j]=max(s[i-1][j]-sigma,s[i][j-1]-sigma,s[i-1][j-1]+score[v[i-1]][w[j-1]])
            if not v[i-1]==w[j-1]:
                s[i][j]=max(s[i-1][j]-sigma,s[i][j-1]-sigma,s[i-1][j-1]-1)
            else:
                s[i][j]=s[i-1][j-1]+score
            if s[i][j]==s[i-1][j]-sigma:
                bk[i][j]=1
            elif s[i][j]==s[i][j-1]-sigma:
                bk[i][j]=2
            else:
                bk[i][j]=3
    AlignV=''
    AlignW=''
    i=len(v)
    j=len(w)
    while True:
        if i>0 and j>0:
            if bk[i][j] ==3:

                AlignV=v[i-1]+AlignV
                AlignW=w[j-1]+AlignW
                i=i-1
                j=j-1
            elif bk[i][j] ==1:
            
                AlignW='-'+AlignW
                AlignV=v[i-1]+AlignV
                i =i-1
            else:
            
                AlignW=w[j-1]+AlignW
                AlignV='-'+AlignV
                j =j-1
        elif i>0:
            AlignW='-'+AlignW
            AlignV=v[i-1]+AlignV
            i =i-1
        elif j>0:
            AlignW=w[j-1]+AlignW
            AlignV='-'+AlignV
            j =j-1
        else:
            break
    mm=''
    m=0
    for i in range(len(AlignV)):
        if not AlignV[i]==AlignW[i]:
            mm=mm+str(i+1)+'('+AlignV[i]+'->'+AlignW[i]+')'+';'
        else:
            m=m+1
    return mm,m/len(w)
def Align(primer,db):

    list_hit=[]




    for seq in SeqIO.parse(db,'fasta'): 
        hit={}     
        hit['stitle']=seq.description
        pos,mm,m=FittingAlignment(str(seq.seq).upper(),primer,1,1)
        hit['qstart']=1
        hit['sstart']=pos-len(primer)+1
        hit['qlen']=len(primer)
        list_hit.append(hit)
    return list_hit

def extensePrimerToHaplotype(ref_genome,primers):
    ref_sequence=''
    primers_extend=[]
    for seq in SeqIO.parse(ref_genome,'fasta'):
        ref_sequence=str(seq.seq).upper()
    for gp in primers:
        #blast Fo:
        ngp=gp
        blast_Fo=blast(gp['primer'][0]['seq'],ref_genome)
        #get best result
        bitscore=0
        for b in blast_Fo:
            if float(b['bitscore'])>bitscore:
                bitscore=float(b['bitscore'])
                #get haplotype by position, extend 20nu
                ht=ref_sequence[int(b['sstart'])-20:int(b['sstart'])+gp['length']+20]
                ngp['haplotype']=ht

        primers_extend.append(ngp)
    return primers_extend

def collectHaplotypeFromCorpus(primers,db):
    ref_sequence={}
    num_seq=0
    for seq in SeqIO.parse(db,'fasta'):
        num_seq=num_seq+1
        ref_sequence[seq.description]=str(seq.seq).upper()
        
    dict_haplotype={}
    dict_domain={}
    check_list={}
    haplotype_set=set()
    name_ht_dict={}
    for gp in primers:
        print('check primer '+gp['name'])
        dict_haplotype[gp['name']]={}
        dict_domain[gp['name']]={}

        blast_haplotype=blast(gp['haplotype'],db)
        check_list[gp['name']]={}
        for b in blast_haplotype:
            ht=''
          
            if b['stitle'] in check_list[gp['name']]:
               if check_list[gp['name']][b['stitle']]>float(b['bitscore']):
                    continue
            position=int(b['sstart'])-1
            if int(b['qstart'])>1:
                position=int(b['sstart'])-int(b['qstart'])
            region=ref_sequence[b['stitle']][position:position+int(b['qlen'])]
            
            dict_domain[gp['name']][b['stitle']]={}
            ht,ps,pe=getAMP(region,gp['primer'][0]['seq'],gp['primer'][2]['seq'])
            
            dict_domain[gp['name']][b['stitle']]['seq']=ht
            
            dict_domain[gp['name']][b['stitle']]['start']=position+ps+1
            dict_domain[gp['name']][b['stitle']]['end']=position+int(b['qlen'])-pe+1
            if ht=='': 
                ht=region
                dict_domain[gp['name']][b['stitle']]['seq']=region
            
                dict_domain[gp['name']][b['stitle']]['start']=position+1
                dict_domain[gp['name']][b['stitle']]['end']=position+int(b['qlen'])+1
            
            
            check_list[gp['name']][b['stitle']]=float(b['bitscore'])
            
            if not ht in haplotype_set:

                newname=gp['name']+'-'+str(len(dict_haplotype[gp['name']].keys())+1)
                dict_haplotype[gp['name']][newname]={}
                dict_haplotype[gp['name']][newname]['seq']=ht
                p,mm,score=FittingAlignment(gp['haplotype'],ht,1,1)
                dict_haplotype[gp['name']][newname]['mm']=str(mm)
                dict_haplotype[gp['name']][newname]['sample']=[]
                haplotype_set.add(ht)
                name_ht_dict[ht]=newname
                dict_haplotype[gp['name']][newname]['sample'].append((b['stitle'],position))
            else:
                dict_haplotype[gp['name']][name_ht_dict[ht]]['sample'].append((b['stitle'],position))
    return dict_haplotype,dict_domain,num_seq
def collectFullDomainFromCorpus(primers,db):
    ref_sequence={}
    for seq in SeqIO.parse(db,'fasta'):
        ref_sequence[seq.description]=str(seq.seq).upper()
    dict_domain={}


    for gp in primers:
        dict_domain[gp['name']]={}

        blast_haplotype_start=blast(gp['primer'][0]['seq'],db)
        blast_haplotype_end=blast(gp['primer'][2]['seq'],db)
        #blast_haplotype_start=Align(gp['primer'][0]['seq'],db)
        #blast_haplotype_end=Align(gp['primer'][2]['seq'],db)

        position_start={}
        position_end={}
        for b in blast_haplotype_start:
            if not b['stitle'] in position_start:
                position_start[ b['stitle']]=[]

            if int(b['qstart'])>1:
                position_start[b['stitle']].append(int(b['sstart'])-int(b['qstart']))
            else:
                position_start[b['stitle']].append(int(b['sstart'])-1)
        for b in blast_haplotype_end:

            if not b['stitle'] in position_end:
                position_end[ b['stitle']]=[]
            position_end[b['stitle']].append(int(b['sstart'])+int(b['qlen'])-1)


        for seq in ref_sequence:
            correct_p=[]
            if seq in position_start and seq in position_end:
                correct_p=getGoodCombination(position_start[seq],position_end[seq],gp['length'])
            dict_domain[gp['name']][seq]={}
            ht=''
            if len(correct_p)>0:
                dict_domain[gp['name']][seq]['seq']=ref_sequence[seq][correct_p[0]:correct_p[1]]
                dict_domain[gp['name']][seq]['start']=correct_p[0]
                dict_domain[gp['name']][seq]['end']=correct_p[1]
            else:
                dict_domain[gp['name']][seq]['seq']=''
                dict_domain[gp['name']][seq]['start']=-1
                dict_domain[gp['name']][seq]['end']=-1



    return dict_domain
def Score(p, q):
    # your code here
    n=len(q)
    count=0
    list_mm=[]
    count_miss=0
    for i in range(0,n):
        if p[i]==q[i]:
            count=count+1
        elif p[i]=='Y' and (q[i]=='C' or q[i]=='T'):
            count=count+1
        elif q[i]=='Y' and (p[i]=='C' or p[i]=='T'):
            count=count+1
        elif p[i]=='R' and (q[i]=='A' or q[i]=='G'):
            count=count+1
        elif q[i]=='R' and (p[i]=='A' or p[i]=='G'):
            count=count+1
        elif p[i]=='S' and (q[i]=='G' or q[i]=='C'):
            count=count+1
        elif q[i]=='S' and (p[i]=='G' or p[i]=='C'):
            count=count+1
        else:
            count_miss=count_miss+1
            list_mm.append(str(i+1)+':'+p[i]+"->"+q[i])
    #if count_miss>3:
     #   list_mm=[]
      #  list_mm.append('not found')
    return count,list_mm
def getMMLong(p, q):
    # your code here
    n=len(q)
    count=0
    list_mm=[]
    count_miss=0
    for i in range(0,n):
        if p[i]==q[i]:
            count=count+1
        else:
            count_miss=count_miss+1
            list_mm.append(str(i+1)+':'+p[i]+"->"+q[i])
    #if count_miss>5:
    #   list_mm=[]
    #  list_mm.append('need alignment')
    return count,list_mm
def getMM(haplotype,primers):
    mm_s=''
    scores=[]
    for primer in primers['primer']:
        k=len(primer['seq'])
        q=str(primer['seq']).upper()
        best=0
        mul_s=''
        list_m=[]
        for i in range(len(haplotype)-k+1):
            p=haplotype[i:i+k]
            score,list_mm=Score(q,p)
            if score>best:
                best=score
                mul_s=''
                for mm in list_mm:
                    mul_s=mul_s+','+mm
        scores.append(best/k)
        if not mul_s=='':
            mm_s=mm_s+primer['type']+'('+mul_s+')'+';'
    return mm_s,scores

def checkAmpliconWithRef(amplicon,ref_region):
    #ret=blast(amplicon,ref_db)
    #get highest score hit
    # bestscore=0
    # best_pos_start=0
    # best_pos_end=0
    # for h in ret:
    #     if float(h['bitscore'])>bestscore:
    #         bestscore=float(h['bitscore'])
    #         best_pos_start=int(h['sstart'])-int(h['qstart'])
    #         best_pos_end=int(h['send'])+(int(h['qlen'])-int(h['qend']))
    
    # ref_amplicon= ref_seq[best_pos_start:best_pos_end]       
    p,mm,m=FittingAlignment(ref_region,amplicon,1,1)
    return mm,m
def export_file(dict_haplotype,primers,output):
    #    blast_fields={'qseqid':t[0], 'qstart':t[1], 'qend':t[2], 'qlen':t[3],\
    #     'sseqid':t[4], 'sstart':t[5], 'send':t[6], 'slen':t[7], 'sstrand':t[8],\
    #      'length':t[10], 'mismatch':t[11], 'qseq':t[12], 'sseq':t[13],\
    #       'stitle':t[14]}
    file_out=os.path.join(output,'haplotype_primer.tsv')
    f=open(file_out,'w')
    f.write('PRIMER\tHAPLOTYPE\tSEQUENCE\tNUMBER\tTOTAL\tFREQUENCY\tMISMATCH WITH PRIMER\tMISMATCH WITH REF\n')
    for gp in primers:
        f.write(gp['name']+'-ref\t\t'+gp['haplotype']+'\t\t\t\t\t\n')
    #statistic_haplotype={}
    statistic_sample={}
    header='SAMPLE ID'
    for gp in dict_haplotype:
        print('export result '+gp)
        total_set=set()
        header=header+'\t'+gp
        for h in dict_haplotype[gp]:
            for sample in dict_haplotype[gp][h]['sample']:
                total_set.add(sample)
                if not sample[0] in statistic_sample:
                    statistic_sample[sample[0]]={}
                statistic_sample[sample[0]][gp]=h
        for h in dict_haplotype[gp]:
            number_haplotype=len(dict_haplotype[gp][h]['sample'])
            total=len(total_set)
            freq=number_haplotype/total
            mm_s=''
            for p in primers:
                if p['name']==gp:
                    mm_s,identity=getMM(dict_haplotype[gp][h]['seq'],p)
            f.write(gp+'\t'+h+'\t'+dict_haplotype[gp][h]['seq']+'\t'+str(number_haplotype)+'\t'+str(total)+'\t'+str(freq)+'\t'+mm_s+'\t'+dict_haplotype[gp][h]['mm']+'\n')
    f.close()
    sample_haplotype_file=os.path.join(output,'sample_haplotype.tsv')
    f=open(sample_haplotype_file,'w')
    f.write(header+'\n')
    set_group_haplotype=set()
    statistic_group_haplotype={}
    for sample in statistic_sample:
        s=sample
        group_haplotype=''
        for gp in dict_haplotype:
            if gp in statistic_sample[sample]:
                s=s+'\t'+statistic_sample[sample][gp]
                group_haplotype=group_haplotype+','+statistic_sample[sample][gp]
            else:
                s=s+'\t'
        if not group_haplotype=='':
            if not group_haplotype[1:] in set_group_haplotype:
                statistic_group_haplotype[group_haplotype[1:]]=[]
                set_group_haplotype.add(group_haplotype[1:])
            statistic_group_haplotype[group_haplotype[1:]].append(sample)
        f.write(s+'\t'+group_haplotype+'\n')
    f.close()
    group_haplotype_file=os.path.join(output,'group_haplotype.tsv')
    f=open(group_haplotype_file,'w')
    f.write('GROUP HAPLOTYPE\tNUMBER\tFREQUENCY\n')
    for group in statistic_group_haplotype:
         f.write(group+'\t'+str(len(statistic_group_haplotype[group]))+'\t'+str(len(statistic_group_haplotype[group])/len(statistic_sample))+'\n')
    f.close()
    haplotype_file=os.path.join(output,'haplotype.fasta')
    f=open(haplotype_file,'w')
    for gp in dict_haplotype:
        for h in dict_haplotype[gp]:
            f.write('>'+h+'\n')
            f.write(dict_haplotype[gp][h]['seq']+'\n')
    f.close()
def export_domain_file(dict_domain,primers,db,file_out,ref_db):
    #    blast_fields={'qseqid':t[0], 'qstart':t[1], 'qend':t[2], 'qlen':t[3],\
    #     'sseqid':t[4], 'sstart':t[5], 'send':t[6], 'slen':t[7], 'sstrand':t[8],\
    #      'length':t[10], 'mismatch':t[11], 'qseq':t[12], 'sseq':t[13],\
    #       'stitle':t[14]}
    #just get 1 seq
    ref_seq=''
    for seq in SeqIO.parse(ref_db,'fasta'):
        ref_seq=str(seq.seq).upper()
        break
    f=open(file_out,'w')

    f.write('SAMPLE\tGROUP PRIMERS\tSTART\tEND\tREGION\tMISMATCH\tIDENTITY\t\t\t\tAMPLICON\tIDENT-REF\tMM REF\tIS MISS')
    f.write('\n')
    set_missAll=set()
    count_sample=0
    for seq in SeqIO.parse(db,'fasta'):
        count_sample=count_sample+1
        isMissAll=True
        for gp in primers:
           
           
            if seq.description in dict_domain[gp['name']] and not dict_domain[gp['name']][seq.description]['seq']=='':
                mm_s,identity=getMM(dict_domain[gp['name']][seq.description]['seq'],gp)
                isMiss=''
                idens=''
                for ide in identity:
                    if ide<0.8:
                        isMiss='Miss'
                    
                    idens=idens+str(ide)+'\t'
                amplicon,ps,pe=getAMP(dict_domain[gp['name']][seq.description]['seq'],gp['primer'][1]['seq'],gp['primer'][2]['seq'])
                if not amplicon=='':
                    mm_r,mr=checkAmpliconWithRef(amplicon,gp['haplotype'])
                    f.write(seq.description+'\t'+gp['name']+\
                    '\t'+str(dict_domain[gp['name']][seq.description]['start'])+\
                    '\t'+str(dict_domain[gp['name']][seq.description]['end'])+\
                    '\t'+dict_domain[gp['name']][seq.description]['seq']+\
                        '\t'+mm_s+'\t'+str(idens)+'\t'+amplicon+'\t'+str(mr)+'\t'+mm_r+'\t'+isMiss+'\n')
                else:
                    isMiss='Miss'
                    f.write(seq.description+'\t'+gp['name']+\
                    '\t'+str(dict_domain[gp['name']][seq.description]['start'])+\
                    '\t'+str(dict_domain[gp['name']][seq.description]['end'])+\
                    '\t'+dict_domain[gp['name']][seq.description]['seq']+\
                        '\t'+mm_s+'\t'+str(idens)+'\t'+amplicon+'\t\t\tMiss\n')
                    
            else:
                 isMiss='Miss'
                 f.write(seq.description+'\t'+gp['name']+\
                '\t'+(str(-1))+\
                '\t'+(str(-1))+\
                    '\t'+\
                        '\t\t\t\t\t\tMiss\n')
            if not isMiss=='Miss':
                isMissAll=False
        if isMissAll==True:
           set_missAll.add(seq.description) 
    f.close()
    # f=open('/media/ktht/Store/Quang/bio/sample_miss_all_primer'+'.tsv','w')
    # for sample in set_missAll :
    #    f.write(sample+'/n')
    # f.close()
def readPrimerFile(primer_file):
    primers=[]
    file1 = open(primer_file, 'r')
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
            p['type']=tok[1]
            p['seq']=tok[2].strip()
            if p['order']<0:
                p['seq']=ReverseComplement(p['seq'])
            primer['primer'].append(p)
        primers.append(primer)
    return primers

   
# def checkAmpliconWithRef(amplicon,ref_db, ref_seq):
#     ret=blast(amplicon,ref_db)
#     #get highest score hit
#     bestscore=0
#     best_pos_start=0
#     best_pos_end=0
#     for h in ret:
#         if float(h['bitscore'])>bestscore:
#             bestscore=float(h['bitscore'])
#             best_pos_start=int(h['sstart'])-int(h['qstart'])
#             best_pos_end=int(h['send'])+(int(h['qlen'])-int(h['qend']))
    
#     ref_amplicon= ref_seq[best_pos_start:best_pos_end]       
#     mm,m=GlobalAlignment(amplicon,ref_amplicon,1,5)
#     return mm,m
def readHaplotypeFile(haplotype_file):
    dict_domain={}
    with open(haplotype_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')

        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                print(f'Column names are {", ".join(row)}')
                line_count += 1
            else:
                if  not row[1] in dict_domain:
                    dict_domain[row[1]]={}
                if not row[0] in dict_domain[row[1]]:
                    dict_domain[row[1]][row[0]]={}
                dict_domain[row[1]][row[0]]['seq']=row[4]
                dict_domain[row[1]][row[0]]['start']=row[2]
                dict_domain[row[1]][row[0]]['end']=row[3]
                line_count += 1
    return dict_domain
def multithread(num,dict_domain,primers,db_file,output,ref):
    print('start thread '+str(num))
    outfile=os.path.join(output,'domain_primer_'+str(num)+'.tsv')
    export_domain_file(dict_domain,primers,db_file,outfile,ref)
    print('end thread '+str(num))
from multiprocessing import Process
def pipeline(args):
    db_file=args.db
    primer_file=args.primer
    primers=readPrimerFile(primer_file)
    setupdbRef(args.ref)
    primers_extend=extensePrimerToHaplotype(args.ref,primers)
    dict_haplotype,dict_domain,num_seq=collectHaplotypeFromCorpus(primers_extend,db_file)
    with open(f"{args.output}/dump_domain.json", 'w') as outfile:
        json.dump(dict_domain, outfile)
    with open(f"{args.output}/dump_dict_haplotype.json", 'w') as outfile:
        json.dump(dict_haplotype, outfile)
    with open(f"{args.output}/dump_num_seq.json", 'w') as outfile:
        json.dump(num_seq, outfile)

    # export_file(dict_haplotype,primers,args.output)
    # dict_haplotype={}
    # dict_domain={}
    # num_seq=0
    with open(f"{args.output}/dump_num_seq.json") as infile:
        num_seq=json.load( infile) 
    with open(f"{args.output}/dump_domain.json") as infile:
        dict_domain=json.load( infile) 
    with open(f"{args.output}/dump_dict_haplotype.json") as infile:
        dict_haplotype=json.load( infile)    
    #export_domain_file(dict_domain,primers,db_file,'/media/ktht/Store/Quang/bio/domain_primer_ultramp_gisaid70k_1N.tsv','/media/ktht/Store/Quang/bio/MN908947.fasta')
    #split db
    # sp_size=10000
    count=0
    num=0
    # handle = None
    batch_size=num_seq/args.threads
    for seq in SeqIO.parse(db_file,'fasta'):
        if count==0:
            num=num+1
            sequence_file_cut=os.path.join(args.output,'seq_'+str(num)+'.fasta')
            handle = open(sequence_file_cut,"w")
        SeqIO.write(seq,handle,"fasta")
        count=count+1
        if count>batch_size:
            count=0
    threads = args.threads
    procs = []
    for n in range(threads):
        if(n>=6):
            sequence_file_cut=os.path.join(args.output,'seq_'+str(n+1)+'.fasta')
            db_file_temp=setupdb(sequence_file_cut)
            proc = Process(target=multithread,args=(n+1,dict_domain,primers,db_file_temp,args.output,args.ref))
            procs.append(proc)
            proc.start()
    for proc in procs:
        proc.join()
    
def main(arguments=sys.argv[1:]):
    # #read primer text:
    # db_file='/media/ktht/Store/Quang/bio/gisaid70k.fasta'
    # db_file_filtered='/media/ktht/Store/Quang/bio/gisaid70k_1N.fasta'
    # db_file=filterFastaFile(db_file,db_file_filtered)
    # db_file=setupdb(db_file)
    # primers=readPrimerFile('/media/ktht/Store/Quang/bio/primerUltramp.txt')
    # setupdbRef('/media/ktht/Store/Quang/bio/MN908947.fasta')
    # primers_extend=extensePrimerToHaplotype('/media/ktht/Store/Quang/bio/MN908947.fasta',primers)
    # dict_haplotype,dict_domain=collectHaplotypeFromCorpus(primers_extend,db_file)
    # export_file(dict_haplotype,primers,'/media/ktht/Store/Quang/bio/haplotype_primer_Ultramp_gisaid70k_1N.tsv')
    # export_domain_file(dict_domain,primers,db_file,'/media/ktht/Store/Quang/bio/domain_primer_ultramp_gisaid70k_1N.tsv','/media/ktht/Store/Quang/bio/MN908947.fasta')
    # db_file='/media/ktht/Store/Quang/bio/gisaid70k.fasta'
    # db_file_filtered='/media/ktht/Store/Quang/bio/gisaid70k_1N.fasta'
    # db_file=filterFastaFile(db_file,db_file_filtered)
    # db_file=setupdb(db_file)
    # primers=readPrimerFile('/media/ktht/Store/Quang/bio/primerUltramp.txt')
    # setupdbRef('/media/ktht/Store/Quang/bio/MN908947.fasta')
    # primers_extend=extensePrimerToHaplotype('/media/ktht/Store/Quang/bio/MN908947.fasta',primers)
    # dict_haplotype,dict_domain=collectHaplotypeFromCorpus(primers_extend,db_file)
    # export_file(dict_haplotype,primers,'/media/ktht/Store/Quang/bio/haplotype_primer_Ultramp_gisaid70k_1N.tsv')
    # #export_domain_file(dict_domain,primers,db_file,'/media/ktht/Store/Quang/bio/domain_primer_ultramp_gisaid70k_1N.tsv','/media/ktht/Store/Quang/bio/MN908947.fasta')
    # #split db
    # # sp_size=10000
    # count=0
    # num=0
    # # handle = None
    # for seq in SeqIO.parse(db_file,'fasta'):
    #     if count==0:
    #         num=num+1
    #         handle = open('/media/ktht/Store/Quang/bio/sequences'+str(num)+".fasta","w")
    #     SeqIO.write(seq,handle,"fasta")
    #     count=count+1
    #     if count==7000:
    #         count=0
    # numbs = [1,2,3,4,5,6,7,8,9,10,11]
    # procs = []
    # for n in range(num):
    #     db_file_temp=setupdb('/media/ktht/Store/Quang/bio/sequences'+str(n+1)+'.fasta')
    #     proc = Process(target=multithread,args=(n+1,dict_domain,primers,db_file_temp))
    #     procs.append(proc)
    #     proc.start()
    # for proc in procs:
    #     proc.join()

    
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
    run_cmd.add_argument('--ref', help='reference sample', type=str) 
    run_cmd.add_argument('--threads', help='Number of thread', type=int) 

    args = parser.parse_args(arguments)
    return args.func(args)
if __name__ == "__main__":
    main()
