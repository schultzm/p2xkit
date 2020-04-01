# bioemboss

Uses Biopython emboss primerscan to align primers to template to find amplimers, and aligns qPCR probes to amplimers.  
Uses bowtie2 to map the probes  
Uses Biopython to form the alignment  
Plus some custom scripting  

To do:  

- develop into a standalone package  
- Refactor code  
- CI testing  
- Testsuite  
- parse args with ArgParse instead of sys.argv  

To use:  
`chmod +x bioemboss.py`  
install dependencies (todo, setup.py)  
`./bioemboss/bioemboss.py test.fa primers_IUPAC.tab probes.fa testoutdir`  

In the command above:  

- `test.fa` is the reference template strand.  
- `primers_IUPAC.tab` is the tab-delimited primer-table as needed by primerscan (EMBOSS)  
	- first column is primer-pair name  
	- second column is forward primer sequence  
	- third column is reverse primer sequence  
	- no header row  
- `probes.fa` is the fasta file of probe sequences  
- `testoutdir` is the parent outdir (must exist prior to running)  

Output is a directory tree containing:  

```

testoutdir/template_name
├── primer_01
│   └── probe_0
│       └── amplimer_0.fa
├── primer_02
│   └── probe_0
│       └── amplimer_0.fa
├── primer_03
│   └── probe_0
│       └── amplimer_0.fa
├── primer_04
│   └── probe_0
│       └── amplimer_0.fa
├── primer_05
│   └── probe_0
│       └── amplimer_0.fa
├── primer_06
│   └── probe_0
│       └── amplimer_0.fa
├── primer_07
│   └── probe_0
│       └── amplimer_0.fa
│       └── amplimer_1.fa
│   └── probe_1
│       └── amplimer_0.fa
├── primer_08
│   └── probe_0
│       └── amplimer_0.fa

```

Result in a fasta file, for primer_11, probe_0, Amplimer 0, looks like:  

```
>templateid | primer_11~~~Amplimer 0~~~28610:28709
GGGGAACTTCTCCTGCTAGAATGGCTGGCAATGGCGGTGATGCTGCTCTTGCTTTGCTGC
TGCTTGACAGATTGAACCAGCTTGAGAGCAAAATGTCTG
>forward_oligo_primer_11 primer_11
GGGGAACTTCTCCTGCTAGAAT--------------------------------------
---------------------------------------
>reverse_oligo_primer_11 primer_11
------------------------------------------------------------
-----------------CAGCTTGAGAGCAAAATGTCTG
>probe_0_primer_11 primer_11
-----------------------------------------------------TTGCTGC
TGCTTGACAGATT--------------------------

```

The coordinates of the hit are given in the templateid fasta header.  