import subprocess
import os, shutil, glob
import re
import sys
import itertools
from Bio import SeqIO
def calVariantsBwa(ref_fa,query_fa,output):
    '''Index reference genome, call aligment with BWA and convert to vcf'''
    temp_folder=output
    
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
    cmd = 'samtools view -hF 2308  {output}/aln-se.sam > {output}/aln-filtered.sam'.format(
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
def parseVariantsFileToList(variant_file,basename):
    f = open(variant_file)
    # use readline() to read the first line
    list_variants=[]
    line = f.readline()
    
    while line:

        if not line.startswith('#'):
            print(line)
            token=line.split('\t')
            obj={'Sample':basename,'GeneID':token[0],'POS':token[1],'REF':token[3],'ALT':token[4]}
            list_variants.append(obj)
        line = f.readline()
    f.close()
    return list_variants
def analysisVariants(variant_file):
    f=open(variant_file,'r')
    line = f.readline()
    list_variants={}
    dict_mutation={}
    while line:

        print(line)
        token=line.split('\t')
        obj={'Sample':token[0],'POS':token[1],'REF':token[2],'ALT':token[3].strip(),'Mutation':token[1].strip()+':'+token[2].strip()+'-'+token[3].strip()}
        if not obj['Mutation'].strip() in dict_mutation:
            dict_mutation[obj['Mutation'].strip()]=set()
        dict_mutation[obj['Mutation'].strip()].add(obj['Sample'])
        if not obj['Sample'] in list_variants:
            list_variants[obj['Sample']]=[]
        list_variants[obj['Sample']].append(obj['Mutation'])
   
        line = f.readline()
    f.close()
    unique_mutation=set()
    f=open('/mnt/data/coronacheck/unique_mutation.tsv','w')
    f.write('Mutation\tCount\n')
    for mu in dict_mutation:
        f.write(mu+'\t'+str(len(dict_mutation[mu])) +'\n')
        if len(dict_mutation[mu])==1:
            unique_mutation.add(mu)
    f.close()
    f=open('/mnt/data/coronacheck/sample_unique_mutation.tsv','w')
    f.write('Sample\tUnique Mutation\tPercent\n')
    for sample in list_variants:
        count_unique=0
        for variant in list_variants[sample]:
            if variant in unique_mutation:
                count_unique=count_unique+1
        perc=str(count_unique/len(dict_mutation)*100)+'%'
        f.write(sample+'\t'+str(count_unique)+'\t'+perc+'\n')
    f.close()
    


def main():
    # ref='/home/quang/Downloads/gisaid4000/EPI_ISL_402124.fasta'
    # queries='/mnt/data/coronacheck/sarscov2_gisaid6000_1N.fasta'
   
    # list_variants=[]
    # for seq in SeqIO.parse(queries,'fasta'): 
    #     tempfolder='/mnt/data/coronacheck/variants/temp'
    #     if not os.path.exists(tempfolder):
    #         os.makedirs(tempfolder)
    #     SeqIO.write(seq,tempfolder+'/query.fasta','fasta')
    #     calVariantsBwa(ref,tempfolder+'/query.fasta',tempfolder)
    #     variants=parseVariantsFileToList(tempfolder+'/variants.vcf',seq.description)
    #     list_variants.extend(variants)
    #     if os.path.exists(tempfolder):
    #         shutil.rmtree(tempfolder)
    # f=open('/mnt/data/coronacheck/variants.tsv','w')
    # f.write('Sample\tPos\tREF\tALT\n')
    # for v in list_variants:
    #     f.write(v['Sample']+'\t'+str(v['POS'])+'\t'+v['REF']+'\t'+v['ALT']+'\n')
    # f.close()
    analysisVariants('/mnt/data/coronacheck/variants.tsv')
if __name__ == '__main__':
    main()
