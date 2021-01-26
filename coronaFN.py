import subprocess
import os, shutil, glob
import re
import sys
import argparse
from Bio import SeqIO
def count_gap(text):
    if text.strip() == "": # To take of care of all space input
        return 0
    #count = sum([1 if char=='N' else 0 for char in text ])
    count = text.count('N') # Your error is here, Only check for 1 space instead of 3 spaces
    count =count+ text.count('n')
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
def setupdb(args):
    """
    make blast database from fasta file in db folder,
    :param : fasta file (with folder'folder is the name of db and filename is 'sequences')
    :return:
    """
    #seqfile='/mnt/data/coronacheck/sarscov2.fasta'
    #gisaid_dir='/home/quang/Downloads/gisaid4000'
    input_dir=args.samples
    output_dir=args.output
    threadshol_n=args.threadshold
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print('Setup database...')
    list_seq=[]
   
    for root, dirs, files in os.walk(input_dir):
        for _file in files:
            if _file.endswith(('.fasta')):
                #print(str(root)+'/'+_file)
                for seq in SeqIO.parse(str(root)+'/'+_file,'fasta'):
                    perc_gap=count_gap(seq.seq)  
                    if perc_gap<=threadshol_n and len(seq.seq)>20000:
                       list_seq.append(seq)
    outfile=os.path.join(output_dir,'sequences_'+str(threadshol_n)+'N.fasta')
    SeqIO.write(list_seq,outfile,'fasta')
    print(len(list_seq))
    cmd="makeblastdb -in {path} -title {name} -dbtype {type} -logfile /dev/null".format(
                            path=outfile,
                            name='corona',
                            type='nucl'

    )
    print (cmd)
    os.system(cmd)
def setupdbfile(dbfile):
    """
    make blast database from fasta file in db folder,
    :param : fasta file (with folder'folder is the name of db and filename is 'sequences')
    :return:
    """

    cmd="makeblastdb -in {path} -title {name} -dbtype {type} -logfile /dev/null".format(
                            path=dbfile,
                            name='corona',
                            type='nucl'

    )
    print (cmd)
    os.system(cmd)
def export_file(sample,db,result,output,dict_cds):
    #    blast_fields={'qseqid':t[0], 'qstart':t[1], 'qend':t[2], 'qlen':t[3],\
    #     'sseqid':t[4], 'sstart':t[5], 'send':t[6], 'slen':t[7], 'sstrand':t[8],\
    #      'length':t[10], 'mismatch':t[11], 'qseq':t[12], 'sseq':t[13],\
    #       'stitle':t[14]}
    f=open(output,'w')
    f.write('PRIMER\tPRIMER START\tPRIMER END\tPRIMER LEN\tSAMPLE ID\tSAMPLE START\tSAMPLE END\tSAMPLE LENGTH\tSAMPLE STRAND\
    \tALIGNMENT LENGTH\t%COVER\tIDENTITY\tMISMATCH\tMISMATCH POSITION\tPRIMER ALIGNED SEQUENCE\tSAMPLE SEQUENCE\tSAMPLE TITLE\tGENE\n')
    #print(result)
    dict_sample_primer={}
    for s in result:
        #filter missmatch and not cover:
        isPrint=False
        if not s['sseqid'] in dict_sample_primer:
            dict_sample_primer[s['sseqid']]={}
        mismatch_c=''
        if not s['qseqid'] in dict_sample_primer[s['sseqid']]:
            dict_sample_primer[s['sseqid']][s['qseqid']]={}
            dict_sample_primer[s['sseqid']][s['qseqid']]['c']=''
            dict_sample_primer[s['sseqid']][s['qseqid']]['mm']=99999
            dict_sample_primer[s['sseqid']][s['qseqid']]['range']=''
        if int(s['mismatch'])>0:
            isPrint=True
            for i in range(len(s['sseq'])):
                if not s['qseq'][i]==s['sseq'][i]:
                    mismatch_c=mismatch_c+str(i+1)+':'+s['qseq'][i]+'-'+s['sseq'][i]+','
                   
            dict_sample_primer[s['sseqid']][s['qseqid']]['mm']= int(s['mismatch'])           
            dict_sample_primer[s['sseqid']][s['qseqid']]['c']='('+mismatch_c+')'

        if float(s['pident'])<100:
            isPrint=True
        if (int(s['qlen'])-int(s['length']))>3 :
             isPrint=False
        elif(int(s['qlen'])-int(s['length']))>0:
            isPrint=True
        else:
            if int(s['mismatch'])==0:
                isPrint=False
        #isPrint=True
      
       
            

        if (int(s['qlen'])-int(s['length'])) < dict_sample_primer[s['sseqid']][s['qseqid']]['mm']:
            if dict_sample_primer[s['sseqid']][s['qseqid']]['c']=='':
                dict_sample_primer[s['sseqid']][s['qseqid']]['mm']=int(s['qlen'])-int(s['length'])
            if (int(s['qlen'])-int(s['length']))>0:
                dict_sample_primer[s['sseqid']][s['qseqid']]['range']='('+str(s['qstart'])+'->'+str(s['qend'])+'/'+s['qlen']+','+str(s['sstart'])+'->'+str(s['send'])+')'
            
        
        if isPrint:
            #find genes:
            gene=''
            if not dict_cds==None:
               
                for key in dict_cds:
                    if s['sseqid'] in key:
                        for cds in dict_cds[key]:
                            if int(s['sstart']) >= dict_cds[key][cds]['cds_f'] and int(s['sstart']) <= dict_cds[key][cds]['cds_e']:
                                gene= dict_cds[key][cds]['gene']
            f.write(s['qseqid']+'\t'+str(s['qstart'])+'\t'+str(s['qend'])+'\t'+str(s['qlen'])+'\t'+\
            s['sseqid']+'\t'+str(s['sstart'])+'\t'+str(s['send'])+'\t'+str(s['slen'])+'\t'+\
            str(s['sstrand'])+'\t'+str(s['length'])+'\t'+str((int(int(s['length'])/int(s['qlen'])*100)))+'\t'+str(s['pident'])+'\t'+str(s['mismatch'])+'\t'+mismatch_c+'\t'+s['qseq']+\
            '\t'+s['sseq']+'\t'+s['stitle']+'\t'+gene+'\n')
    f.close()
    return dict_sample_primer

def blast(sample,db, identity=90, threads=1, mincov=90,dbtype='nucl'):
    """
    Call blastn with params
    :param query_file (in fasta), db (blast indexed db), number of threads and identity
    :return: list BLASTFields objects
    """
    #check db is indexed
    #dbfile=os.path.join(db_folder, 'sequences')




    # run blastn
    cmd ='blastn -query {query} -task blastn-short -max_target_seqs 1000000 -perc_identity {identity} -db {db} -outfmt \'6 qseqid qstart qend qlen sseqid sstart send slen sstrand length pident mismatch qseq sseq stitle\'  -num_threads {threads} >temp.tab'.format(

        query=sample,
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
           'stitle':t[14]}
        result.append(blast_fields)
        line = f.readline()
    f.close()
    if os.path.exists('temp.tab'):
        os.remove('temp.tab')



    return result


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



def calVariantsBwa(ref_fa,query_fa,output):
    '''Index reference genome, call aligment with BWA and convert to vcf'''
    temp_folder=output+"/temp"
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)
    #indexing reference
    isIndex=False
    for root, dirs, files in os.walk(temp_folder):
        for _file in files:
            if _file.endswith(('.bwt')):
                isIndex=True
                break

    if not isIndex:
        cmd='bwa index '+ref_fa
        os.system(cmd)
        cmd='samtools faidx '+ref_fa
        os.system(cmd)
    #alignment
    cmd = 'bwa mem {ref} {query} > {output}/aln-se.sam'.format(
        ref=ref_fa,
        query=query_fa,
        output=temp_folder
    )
    os.system(cmd)
    #filter duplicated aligment
    cmd = 'samtools view -hF 2308 {output}/aln-se.sam > {output}/aln-filtered.sam'.format(
        output=temp_folder
    )
    os.system(cmd)
    #convert to bam
    cmd = 'samtools view -S -b {output}/aln-filtered.sam > {output}/aln-filtered.bam'.format(
        output=temp_folder
    )
    os.system(cmd)
    #sort BAM for SNP calling
    cmd = 'samtools sort {output}/aln-filtered.bam > {output}/aln-filtered-sorted.bam'.format(
        output=temp_folder
    )
    os.system(cmd)
    #generate raw bcf
    cmd = 'samtools mpileup -g -f {ref} {output}/aln-filtered-sorted.bam | bcftools call -mv -Ov > {output}/variants.vcf'.format(
        ref=ref_fa,
        output=temp_folder
    )
    os.system(cmd)
    #generate vcf
    '''cmd = 'bcftools view -bvcg {output}/raw.bcf > {output}/var.bcf'.format(

        output=temp_folder
    )
    os.system(cmd)
    #generate vcf
    cmd = 'bcftools view {output}/var.bcf | vcfutils.pl varFilter - > {output}/var-final.vcf'.format(

            output=temp_folder
        )
    os.system(cmd)'''
    return  temp_folder+'/variants.vcf'
def readCDSFile(cds_file):
    dic_cds={}
    for seq in SeqIO.parse(cds_file,'fasta'):
        if seq.description.startswith('join'):
            cdses=seq.description.split('|')
            for i  in range(0,len(cdses)-1):
                cds=cdses[i].replace('join','').replace('(','').replace(')','')
                if not ':' in cds:
                    continue
                id=cds.split(':')[0].strip()
                c=cds.split(':')[1].strip()
                gene=cdses[len(cdses)-1]
                if not id in dic_cds:
                    dic_cds[id]={} 
                dic_cds[id][c]={}
                dic_cds[id][c]['cds_f']=int(c.split('..')[0])
                dic_cds[id][c]['cds_e']=int(c.split('..')[1])
                dic_cds[id][c]['gene']=gene
        else:
            id=seq.name.split(':')[0].strip()
            cds=seq.name.split(':')[1].strip()
            gene=seq.description.split('|')[1]
            if not id in dic_cds:
                dic_cds[id]={} 
            dic_cds[id][cds]={}
        
            dic_cds[id][cds]['cds_f']=int(cds.split('..')[0])
            dic_cds[id][cds]['cds_e']=int(cds.split('..')[1])
            dic_cds[id][cds]['gene']=gene

        
    return dic_cds
def pipeline(args):
    primer_file=args.primer
    db=args.db
    output=args.output
    primers=[]
    file1 = open(primer_file, 'r')
    num_primer = int(file1.readline())
    if not os.path.exists(output):
        os.makedirs(output)
    print('read primer file')

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
    primer_fasta=os.path.join(output,'primer.fasta')
    f=open(primer_fasta,'w')
    for p in primers:
        for primer in p['primer']:
            f.write('>'+p['name']+"|"+str(primer['order'])+"\n")
            f.write(primer['seq']+'\n')
    f.close()
    #setupdb()
    #setupdbfile('/media/ktht/Store/Quang/bio/sequences.fasta')
    #dict_cds=readCDSFile('/media/ktht/Store/Quang/bio/sequences.fasta')
    print('blast primer with db')
    ret=blast(primer_fasta,db=db,mincov=70, identity=70,  threads=8)
    #print(ret)
    outfile=os.path.join(output,'blasthit_primer_FN.tsv')
    print('make sumary file')
    dict_sample_primer=export_file(primer_fasta,'corona',ret,outfile,None)
    outfile_summary=os.path.join(output,'summary.tsv')
    f=open(outfile_summary,'w')
    header='Sample'
    for p in primers:
            for primer in p['primer']:
                header=header+'\t'+p['name']+'|'+str(primer['order'])
    f.write(header+'\n')
    for sample in dict_sample_primer:
        s=sample
        for p in primers:
            for primer in p['primer']:
                pt=p['name']+'|'+str(primer['order'])
                print(pt)
                print(dict_sample_primer[sample])
                if pt in dict_sample_primer[sample]:
                    s=s+'\t'+str(dict_sample_primer[sample][pt]['mm'])+dict_sample_primer[sample][pt]['c']+dict_sample_primer[sample][pt]['range']
                else:
                    s=s+'\t'+'not found'
        f.write(s+'\n')
    f.close()
def main(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        prog='corona_fn_checker',
        description='Tool for checking False Negative error when testing specified primers with corona')
    subparsers = parser.add_subparsers(title='sub command', help='sub command help')    
    setup_cmd = subparsers.add_parser(
        'setupdb', description='Collect and filter from corona fasta files', help='Setup corona db file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    setup_cmd.set_defaults(func=setupdb)   
    setup_cmd.add_argument('--samples', help='Directory of collection of samples in fasta format',type=str)
    setup_cmd.add_argument('--output', help='Directory of corona db file',type=str)
    setup_cmd.add_argument('--threadshold', help='Percentage of N character threadshold in sample', type=int) 
    run_cmd = subparsers.add_parser(
        'run', description='Check false postive FN with specified primer', help='Check FN',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    run_cmd.set_defaults(func=pipeline)   
    run_cmd.add_argument('--primer', help='Primer list in text file',type=str)
    run_cmd.add_argument('--db', help='Corona db file',type=str)
    run_cmd.add_argument('--output', help='Output folder', type=str) 
    args = parser.parse_args(arguments)
    return args.func(args)
if __name__ == "__main__":
    main()
