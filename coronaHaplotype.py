import subprocess
import os, shutil, glob
import re
import sys
import itertools
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
    
    # list_seq=[]
   
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
def blast(sample,db, identity=90, threads=1, mincov=90,dbtype='nucl'):
    """
    Call blastn with params
    :param query_file (in fasta), db (blast indexed db), number of threads and identity
    :return: list BLASTFields objects
    """
    #check db is indexed
    #dbfile=os.path.join(db_folder, 'sequences')
    # run blastn
    f=open('temp.fasta','w')
    f.write('>primer\n')
    f.write(sample)
    f.close()
    cmd ='blastn -query {query} -task blastn-short -max_target_seqs 1000000 -perc_identity {identity} -db {db} -outfmt \'6 qseqid qstart qend qlen sseqid sstart send slen sstrand length pident mismatch qseq sseq stitle bitscore\'  -num_threads {threads} >temp.tab'.format(

        query='temp.fasta',
        identity=identity,
        db=db,
        threads=threads

    )
    print(cmd)
    os.system(cmd)
    #parse result
    f=open('temp.tab')
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
    if os.path.exists('temp.tab'):
        os.remove('temp.tab')
    if os.path.exists('temp.fasta'):
        os.remove('temp.fasta')
    return result

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
    for seq in SeqIO.parse(db,'fasta'):
        ref_sequence[seq.description]=str(seq.seq).upper()
    dict_haplotype={}
    check_list={}
    haplotype_set=set()
    name_ht_dict={}
    for gp in primers:
        dict_haplotype[gp['name']]={}
        
        blast_haplotype=blast(gp['haplotype'],db)
        check_list[gp['name']]={}
        for b in blast_haplotype:
            ht=''
            isNeedUpdate=False
            if b['stitle'] in check_list[gp['name']]:
               if check_list[gp['name']][b['stitle']]>float(b['bitscore']):
                    continue
            position=int(b['sstart'])-1
            if int(b['qstart'])>1:
                position=int(b['sstart'])-int(b['qstart'])
            ht=ref_sequence[b['stitle']][position:position+int(b['qlen'])]
            check_list[gp['name']][b['stitle']]=float(b['bitscore'])
            if not ht in haplotype_set:
            
                newname=gp['name']+'-'+str(len(dict_haplotype[gp['name']].keys())+1)
                dict_haplotype[gp['name']][newname]={}
                dict_haplotype[gp['name']][newname]['seq']=ht
                score,mm=getMMLong(gp['haplotype'],ht)
                dict_haplotype[gp['name']][newname]['mm']=str(mm)
                dict_haplotype[gp['name']][newname]['sample']=[]
                haplotype_set.add(ht)
                name_ht_dict[ht]=newname
                dict_haplotype[gp['name']][newname]['sample'].append((b['stitle'],position))
            else:
                dict_haplotype[gp['name']][name_ht_dict[ht]]['sample'].append((b['stitle'],position))
    return dict_haplotype
def collectFullDomainFromCorpus(primers,db):
    ref_sequence={}
    for seq in SeqIO.parse(db,'fasta'):
        ref_sequence[seq.description]=str(seq.seq).upper()
    dict_domain={}

   
    for gp in primers:
        dict_domain[gp['name']]={}
        
        blast_haplotype_start=blast(gp['primer'][0]['seq'],db)
        blast_haplotype_end=blast(gp['primer'][2]['seq'],db)
        
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
            position_end[b['stitle']].append(int(b['sstart'])+int(b['qlen']))         
    

        for seq in ref_sequence:
            correct_p=getGoodCombination(position_start[seq],position_end[seq],gp['length'])
            dict_domain[gp['name']][seq]={}
            ht=''
            if len(correct_p)>0:
                dict_domain[gp['name']][seq]['seq']=ref_sequence[seq][correct_p[0]:correct_p[1]+1]
                dict_domain[gp['name']][seq]['start']=correct_p[0]
                dict_domain[gp['name']][seq]['end']=correct_p[1]+1
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
    if count_miss>3:
        list_mm=[]
        list_mm.append('not found')
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
    if count_miss>5:
        list_mm=[]
        list_mm.append('need alignment')
    return count,list_mm
def getMM(haplotype,primers):
    mm_s=''
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
        if not mul_s=='':
            mm_s=mm_s+primer['type']+'('+mul_s+')'+';'

    return mm_s

def export_file(dict_haplotype,primers,file_out):
    #    blast_fields={'qseqid':t[0], 'qstart':t[1], 'qend':t[2], 'qlen':t[3],\
    #     'sseqid':t[4], 'sstart':t[5], 'send':t[6], 'slen':t[7], 'sstrand':t[8],\
    #      'length':t[10], 'mismatch':t[11], 'qseq':t[12], 'sseq':t[13],\
    #       'stitle':t[14]}
    f=open(file_out,'w')

    f.write('PRIMER\tHAPLOTYPE\tSEQUENCE\tNUMBER\tTOTAL\tFREQUENCY\tMISMATCH WITH PRIMER\tMISMATCH WITH REF\n')
    for gp in primers:
        f.write(gp['name']+'-ref\t\t'+gp['haplotype']+'\t\t\t\t\t\n')
    statistic_haplotype={}
    statistic_sample={}
    header='SAMPLE ID'
    for gp in dict_haplotype:
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
                    mm_s=getMM(dict_haplotype[gp][h]['seq'],p)
            f.write(gp+'\t'+h+'\t'+dict_haplotype[gp][h]['seq']+'\t'+str(number_haplotype)+'\t'+str(total)+'\t'+str(freq)+'\t'+mm_s+'\t'+dict_haplotype[gp][h]['mm']+'\n')
    f.close()
    f=open('/mnt/data/coronacheck/matrix.tsv','w')
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
    f=open('/mnt/data/coronacheck/group_haplotype.tsv','w')
    f.write('GROUP HAPLOTYPE\tNUMBER\tFREQUENCY\n')  
    for group in statistic_group_haplotype:
         f.write(group+'\t'+str(len(statistic_group_haplotype[group]))+'\t'+str(len(statistic_group_haplotype[group])/len(statistic_sample))+'\n')
    f.close()
    f=open('/mnt/data/coronacheck/haplotype.fasta','w')
    for gp in dict_haplotype:
        for h in dict_haplotype[gp]:
            f.write('>'+h+'\n')
            f.write(dict_haplotype[gp][h]['seq']+'\n')
    f.close()
def export_domain_file(dict_domain,primers,db,file_out):
    #    blast_fields={'qseqid':t[0], 'qstart':t[1], 'qend':t[2], 'qlen':t[3],\
    #     'sseqid':t[4], 'sstart':t[5], 'send':t[6], 'slen':t[7], 'sstrand':t[8],\
    #      'length':t[10], 'mismatch':t[11], 'qseq':t[12], 'sseq':t[13],\
    #       'stitle':t[14]}
    f=open(file_out,'w')

    f.write('SAMPLE\tGROUP PRIMERS\tSTART\tEND\tSEQUENCE\tMISMATCH')
    
    f.write('\n')
    for seq in SeqIO.parse(db,'fasta'):
        
        for gp in primers:
            mm_s=getMM(dict_domain[gp['name']][seq.description],gp)
            f.write(seq.description+'\t'+gp['name']+\
                '\t'+str(dict_domain[gp['name']][seq.description]['start'])+\
                '\t'+str(dict_domain[gp['name']][seq.description]['end'])+\
                    '\t'+dict_domain[gp['name']][seq.description]['seq']+\
                        '\t'+mm_s+'\n')
        
    f.close()
    
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



def main(arguments=sys.argv[1:]):
    #read primer text:
    db_file='/mnt/data/coronacheck/sarscov2_gisaid6000_1N.fasta'
    db_file=setupdb(db_file)
    primers=readPrimerFile('/mnt/data/coronacheck/primerUSA.txt')
    setupdbRef('/home/quang/Downloads/coronavirus/NC_045512.2.fna')
    primers_extend=extensePrimerToHaplotype('/home/quang/Downloads/coronavirus/NC_045512.2.fna',primers)
    dict_haplotype=collectHaplotypeFromCorpus(primers_extend,db_file)
    #print(ret)
    export_file(dict_haplotype,primers,'/mnt/data/coronacheck/haplotype_primer_CDC_gisaid6000_1N.tsv')  
    #calVariantsBwa('/mnt/data/coronacheck/sarscov2.fasta','/mnt/data/coronacheck/primer.fasta','/mnt/data/coronacheck/bwa')
    #dict_domain=collectFullDomainFromCorpus(primers_extend,db_file)
    #export_domain_file(dict_domain,primers,db_file,'/mnt/data/coronacheck/domain_primer_ultramp_gisaid6000_1N.tsv')

if __name__ == "__main__":
    main()
