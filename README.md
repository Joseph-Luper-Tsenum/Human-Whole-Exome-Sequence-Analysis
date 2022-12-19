# Human-Whole-Exome-Sequence-Analysis
Analysis of human whole exome sequence from Finnish female populations was carried out, starting with downloading of raw sequence reads, indexing, alignment of sequence reads to the GRCh19 reference human genome, generating a BAM file, variant calling to generate a VCF file, and variant analysis to identify variants of breast cancer melanoma. The sequence data was gotten from individual exomes from the 1000 Genomes Project publicly available at the 1000 Genomes Data Portal https://www.internationalgenome.org/data-portal/sample  A potentially pathogenic or damaging variant (Gene: AGRN; Name: AGRN Protein; Variant ID: 1-978577-G-C) was selected for further downstream analysis to gain insight into the role of this gene and the consequences of other damaging mutations in this gene as a guide to predict what phenotype this variant may have in a homozygous individual. Finally, the protein structure of the wild-type protein (AGRN Protein) was modeled to identify the location of the variant amino acid on the structural model followed by the functional importance of this variant region of the protein and how this particular amino acid change may affect the protein’s structure or function. 

## Commands used: 

(i) Sequence alignment with bwa mem: bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ERR031897_1.fastq.gz ERR031897_2.fastq.gz > ERR031897.sam


(ii) Cleaning up

```samtools fixmate -O bam ERR031897.sam ERR031897fixmate.bam```


(iii) Sorting

```samtools sort -O bam -o  ERR031897sorted.bam -T /tmp/ERR031897temp ERR031897fixmate.bam```


(iv) Indexing

```samtools index ERR031897sorted.bam```

(v) Variant calling with mpileup

```samtools mpileup -go ERR031897.bcf -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ERR031897sorted.bam```

```bcftools call -vmo z -o ERR031897.vcf.gz  ERR031897.bcf```

(vi) Indexing vcf file

```tabix -p ERR031897.vcf  ERR031897.vcf.gz```

```tabix -p vcf ERR031897.vcf.gz```

(vii) Prepare graphs and statistics to assist in filtering variants

```bcftools stats -F GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -s – ERR031897.vcf.gz > ERR031897.vcf.gz.stats```

(viii) Filtering vcf.gz file

```bcftools filter -O z -o ERR031897.filtered.vcf.gz -s LOWQUAL -i'%QUAL>10' ERR031897.vcf.gz```

## Link to wANNOAR output 

https://u27309186.ct.sendgrid.net/ls/click?upn=umGjyFNc-2F7SdYZClx7o5vANo-2B-2B-2B6TfBp4gzxHGDQAoARn8yJ78smd-2FJ0Aj-2FZkzaxJlnvREf-2FybNhmhmX5ToVaG-2BPgzCTukRkEqtBE6lwuL0-3Dg8S7_GBtR8zpR-2BxraP6R65fiRg55Fu05gPqkKgUjrYLN-2F6qCoThB-2Bk4xCaqUun0PIuyJUn0MwXP-2FVXlUPOc1V2Px1rvikIVYEdkCOgoY4488GmiAQ3RvLJqVkDnY-2FuP7E8M0mQaOA575Lw-2FWvokd3SV44td5oHdjnmFJHK8Z7btfA-2Fbmia026-2F0prejgFEM56-2FNHPQ7uIM2DctlBKlhW62-2BBXeQ-3D-3D

## References

1. Jumper, J et al. Highly accurate protein structure prediction with AlphaFold. Nature (2021).
Varadi, M et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. Nucleic Acids Research (2021).

2. Zhang Y, Dai Y, Han JN, Chen ZH, Ling L, Pu CQ, Cui LY, Huang XS. A Novel AGRN Mutation Leads to Congenital Myasthenic Syndrome Only Affecting Limb-girdle Muscle. Chin Med J (Engl). 2017 Oct 5;130(19):2279-2282. doi: 10.4103/0366-6999.215332. PMID: 28937031; PMCID: PMC5634075.

3. Jacquier, A., Risson, V., Simonet, T. et al. Severe congenital myasthenic syndromes caused by agrin mutations affecting secretion by motoneurons. Acta Neuropathol 144, 707–731 (2022). https://doi.org/10.1007/s00401-022-02475-8

4. Andrew Zhang, Steven Shook. Congenital Myasthenic Syndrome Due to a Novel Compound Heterozygous AGRN Gene Variant (1225). Neurology Apr 2021, 96 (15 Supplement) 1225.

5. Kröger S et al. Agrin in the developing CNS: new roles for a synapse organizer. News Physiol Sci 2002 Oct;17207-212. PMID: 12270958.

6. Fagerberg L, Hallström BM, Oksvold P, Kampf C, Djureinovic D, Odeberg J, Habuka M, Tahmasebpoor S, Danielsson A, Edlund K, Asplund A, Sjöstedt E, Lundberg E, Szigyarto CA, Skogs M, Takanen JO, Berling H, Tegel H, Mulder J, Nilsson P, Schwenk JM, Lindskog C, Danielsson F, Mardinoglu A, Sivertsson A, von Feilitzen K, Forsberg M, Zwahlen M, Olsson I, Navani S, Huss M, Nielsen J, Ponten F, Uhlén M. Analysis of the human tissue-specific expression by genome-wide integration of transcriptomics and antibody-based proteomics. Mol Cell Proteomics. 2014 Feb;13(2):397-406. doi: 10.1074/mcp.M113.035600. Epub 2013 Dec 5. PMID: 24309898; PMCID: PMC3916642.

7. Ohkawara B, Shen X, Selcen D, Nazim M, Bril V, Tarnopolsky MA, Brady L, Fukami S, Amato AA, Yis U, Ohno K, Engel AG. Congenital myasthenic syndrome-associated agrin variants affect clustering of acetylcholine receptors in a domain-specific manner. JCI Insight. 2020 Apr 9;5(7):e132023. doi: 10.1172/jci.insight.132023. PMID: 32271162; PMCID: PMC7205260.


