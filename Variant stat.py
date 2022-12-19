# Cheng

# link: https://samtools.github.io/hts-specs/VCFv4.1.pdf


# ---------------------------------------------------------------
# ref file download: 
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/
# GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz


# ---------------------------------------------------------------
# 1000 genomes exome data download
# https://www.internationalgenome.org/data-portal/sample 

# Filter by technology: exome
# Filter by collection: 1000 genomes phase 3 release
# Filter by population (your choice, but avoid SGDP and HGDP)
# Select a Sample ID 
# Check Sequence for Data types
# Check Exome for Technologies
# Download matching pair – same ID_1 and same ID_2.fastq.gz

# Then I pick this one: (everyone need to pick theirselves)
# ERR016121_1.fastq.gz
# ERR016121_2.fastq.gz
# https://www.internationalgenome.org/data-portal/sample/HG00513


# ---------------------------------------------------------------
# samtools download!!!
# link: https://www.htslib.org/download/
cd samtools-1.x  
./configure --prefix=/Users/chengzhang/Desktop/BIOL6150Project2
make
make install

export PATH=/Users/chengzhang/Desktop/BIOL6150Project2/bin:$PATH    # for sh or bash users


# ---------------------------------------------------------------
# bwa download!!!
git clone https://github.com/lh3/bwa.git
cd bwa; make

export PATH=/Users/chengzhang/Desktop/BIOL6150Project2/bwa:$PATH 


# ---------------------------------------------------------------
# mapping
# link: https://www.htslib.org/workflow/wgs-call.html

# solution now
# TA adivse us to use another ref and index by ourself

# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/
# GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz                              
# 2014-01-11 09:09  833M 

# bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz (no need!!!!)

bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ERR016121_1.fastq ERR016121_2.fastq > lane.sam

head -50 lane.sam




# fna.fai
# fna.gz

# do clean up procedure, sort, index

# NA12345bwamemfixmate.bam -- input file

#Annovar web https://wannovar.wglab.org/



