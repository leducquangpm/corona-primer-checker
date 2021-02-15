# Corona-primer-checker
#### There scripts are used for supporting in analyzing primers with the nCov dataset. 
### How to install:
- Requirement:  
  + Environment: python3.6, Biopython
  + Tools: blastn

## CoronaFN.py
Quick report on presence or absence of primers in nCov samples.
### How to use:
Setup dataset: combine all samples sequences into one fasta file, ignore error or bad sequence, prepare database for blast tool. 
  ```
   python coronaFN.py setupdb --samples <folder of samples in fasta files> --output <output folder> --threadshold <percentage of 'N' character>
   ```
   Example:
   ```
   python coronaFN.py setupdb --samples input/samples --output out --threadshold 1
   ```
  Run pipeline:
   ```
    python coronaFN.py run --primer <primer file> --db <nCOv samples> --output <output folder>
   ```
   - primer file in text format (see 2 example primer files)
   - nCov samples: output of setupdb command above.
   - Output file: 
    + summary.txt: quick report presence and absence of primers on each samples
    + blasthit_primer_FN.tsv: more detail of blast result
  ## CoronaHaplotype.py
Detail analysis of primers in nCov samples and some basic statistic about target domains of primers.
### How to use:
  Run pipeline:
   ```
    python coronaHaplotype.py run --primer <primer file> --db <nCOv samples> --output <output folder> --ref <reference sample> --threads <number of threads>
   ```
   - Reference sample: a ref nCov sequence for checking mismatch of target domain on each samples.
   - Threads: number of CPU core, use for spliting dataset into smaller pieces and run parallel for boosting speed. 
   - Output file: 
    + domain_primer.tsv: extract target region of each primer on each sample and check mutation on target region.
    + haplotype_primer.tsv: list all haplotype on target region
    + sample_haplotype.tsv: mapping sample to haplotype 
