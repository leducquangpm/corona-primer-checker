#!/bin/bash
python coronaFN.py setupdb --samples /media/ktht/Store/Quang/bio/corona-data/GISAID --output /media/ktht/Store/Quang/bio/corona_out_25_1 --threadshold 1
python coronaFN.py run --primer /media/ktht/Store/Quang/bio/primerUltramp.txt --db /media/ktht/Store/Quang/bio/corona_out_25_1/sequences_1N.fasta --output /media/ktht/Store/Quang/bio/corona_out_25_1
