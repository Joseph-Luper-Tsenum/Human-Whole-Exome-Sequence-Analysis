# Human-Whole-Exome-Sequence-Analysis
Analysis of human whole exome sequence from Finnish female populations was carried out, starting with downloading of raw sequence reads, indexing, alignment of sequence reads to the GRCh19 reference human genome, generating a BAM file, variant calling to generate a VCF file, and variant analysis to identify variants of breast cancer melanoma. The sequence data was gotten from individual exomes from the 1000 Genomes Project publicly available at the 1000 Genomes Data Portal https://www.internationalgenome.org/data-portal/sample  A potentially pathogenic or damaging variant (Gene: AGRN; Name: AGRN Protein; Variant ID: 1-978577-G-C) was selected for further downstream analysis to gain insight into the role of this gene and the consequences of other damaging mutations in this gene as a guide to predict what phenotype this variant may have in a homozygous individual. Finally, the protein structure of the wild-type protein (AGRN Protein) was modeled to identify the location of the variant amino acid on the structural model followed by the functional importance of this variant region of the protein and how this particular amino acid change may affect the protein’s structure or function. 

## Variant effect on gene/protein function

AGRN:NM_198576:exon1:c.G127T:p.A43S nonsynonymous is the variant of the gene AGRN and AGRN protein. The role of AGRN in neurons is controlled by proteolytic processing, alternative splicing and glycan binding. AGRN regulates calcium ion homeostasis in neurons and it does this by causing an increase in cytoplasmic calcium ions. AGRN is a component of AGRN-LRP4, a receptor complex that plays an important role in inducing the phosphorylation and activation of MUSK. This activation of MUSK in myotubes produces NMJ. It does this by modulating a number of biological processes such as transcription of specific genes and AchR clustering in the postsynaptic membrane. AGRN variant functions differentially in the central nervous system (CNS). It does this by inhibiting the alpha(3)-subtype of Na+/K+-ATPase and inducing depolarization at CNS synapses. The encoded AGRN protein contains a large amount of epidermal growth factor domains, laminin G, and Kazal type serine protease inhibitorn (Fagerberg et al. 2013). Post-translational modification of AGRN protein results in the addition of glycosaminoglycans and disulfide bonds ((Fagerberg et al. 2013). Zhang and Shook 2021, Zhang et al. 2017 and Jacquier et al. 2022 all reported that AGRN defect causes congenital myasthenic syndromes (CMSs) by inducing the development and maintenance of neuromuscular transmission. Kröger et al 2002 also reported that heparan sulfate proteoglycan agrin plays an important role for the formation, maintenance, and regeneration of the neuromuscular junction. Missense variants such as p.S1180L, p.R1509W, p.G1675S, and p.Y1877D weakens agrin-induced AChR clustering in C2C12 myotubes (Ohkawara et al. 2020).

## Individual Sequence ID and Information

HG00330 data details

Sex: Female

Populations: Finnish in Finland, European Ancestry

Biosample ID: SAME123985

Cell line Source: HG00330 at Coriell

Variant Annotation Tool Used: Wannovar

Sample identifier = Exome Project

File Name = ERR031897.filtered.vcf.gz

File format = vcf4

Reference genome = hg19

Disease model = no filtering

Processed variants = 585537

Phenotype: breast cancer melanoma

Exonic variant list from the wANNOVAR output after filtration, the last filter step with variants left and used = (Total: 5429)

Gene list from the wANNOVAR output, input into Phenolyzer = (Total: 3655)

## Summary table of variants

![image](https://user-images.githubusercontent.com/58364462/208487590-e8c45976-71e6-4577-a5ca-b68b0afa1ee8.png)

## Table of potentially damaging or pathogenic nonsynonymous variants

![image](https://user-images.githubusercontent.com/58364462/208488367-52eb4d38-a244-4a6f-9250-e23e2dac855a.png)

## Allele Frequency in the Population

![image](https://user-images.githubusercontent.com/58364462/208488674-b3a50a59-2812-413f-8efe-f9725f1db65a.png)

## Popmax Allele Frequency

![image](https://user-images.githubusercontent.com/58364462/208488898-f9ed4010-4873-4a03-9c9a-1606f5886b7b.png)

## PAthogenic Variant Information

![image](https://user-images.githubusercontent.com/58364462/208489185-25ff4e2a-b427-4104-bde9-9d981a1a77e2.png)

![image](https://user-images.githubusercontent.com/58364462/208489223-30d7dbd8-09e6-4a65-ba82-08f095ab9c95.png)
Gene and Protein Name, and variant ID (rs number): Gene: AGRN; Name: AGRN Protein; Variant ID: 1-978577-G-C
Variant genotype (DNA change and amino acid change): AGRN:NM_198576:exon1:c.G127T:p.A43S nonsynonymous
Variant frequency in the overall human population: 0.0002642
Population with the highest frequency of the variant: East Asia with allele frequency of 0.002836

## Structural Model of AGRN protein

AlphaFold Protein Structure Database was used in predicting the structure of AGRN protein. The 3D structure model of the protein is shown below

![image](https://user-images.githubusercontent.com/58364462/208490312-2a8ebb76-7088-4c08-ac74-b487c06fccf9.png)

## Visualization of altered amino acid on structural model

![image](https://user-images.githubusercontent.com/58364462/208490559-08dc520f-37fa-4bfc-addb-9ca5501b81e0.png)
There is low confidence interval around the site of the altered amino acid shown in yellow, hence the model is not reliable. This means that the change around the site of the variant amino acid will not affect the function of this protein.


## Commands used: 

(i) Sequence alignment with bwa mem

```bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ERR031897_1.fastq.gz ERR031897_2.fastq.gz > ERR031897.sam```


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


