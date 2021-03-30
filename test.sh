#!/bin/bash
#python coronaFN.py setupdb --samples /media/ktht/Store/Quang/bio/corona-data/GISAID --output /media/ktht/Store/Quang/bio/corona_out_25_1 --threadshold 1
#python coronaFN.py run --primer /media/ktht/Store/Quang/bio/primerUltramp.txt --db /media/ktht/Store/Quang/bio/corona_out_25_1/sequences_1N.fasta --output /media/ktht/Store/Quang/bio/corona_out_25_1
#python coronaHaplotype.py run --primer /media/ktht/Store/Quang/bio/primerUltramp.txt --db /media/ktht/Store/Quang/bio/corona_out_25_1/sequences_1N.fasta --output /media/ktht/Store/Quang/bio/corona_out_25_1 --ref /media/ktht/Store/Quang/bio/MN908947.fasta --thread 10
python coronaBlastSearchPrimer.py setupdb --samples /media/ktht/Store/Quang/bio/corona-data/NCBI --output /media/ktht/Store/Quang/bio/corona_out_12_3 --threadshold 1
python coronaBlastSearchPrimer.py run --primer /media/ktht/Store/Quang/bio/primerUltramp.txt --db /media/ktht/Store/Quang/bio/corona_out_12_3/sequences_1N.fasta --output /media/ktht/Store/Quang/bio/corona_out_12_3
python coronaHaplotype.py run --primer /media/ktht/Store/Quang/bio/primerUltramp.txt --db /media/ktht/Store/Quang/bio/corona_out_12_3/sequences_1N.fasta --output /media/ktht/Store/Quang/bio/corona_out_12_3 --ref /media/ktht/Store/Quang/bio/MN908947.fasta --thread 5
python coronaExactSearchPrimer.py run --primer /media/ktht/Store/Quang/bio/primerUltramp.txt --db /media/ktht/Store/Quang/bio/corona_out_12_3/sequences_1N.fasta --output /media/ktht/Store/Quang/bio/corona_out_12_3
