set -x

# file can be *.gz or *.fastq
fq1="celseq2/simu_celseq2_R1.fq"
fq2="celseq2/simu_celseq2_R2.fq"

# no index1
index1_start=-1 
# trim 2bp of read1
index1_len=2

index2_start=6
index2_len=8
umi_start=0
umi_len=6

# experiment name
expr_name=test_celseq2

out_dir=celseq2

mkdir -p $out_dir
mkdir -p $out_dir/count
mkdir -p $out_dir/stat

unaligned_fq=$out_dir/$expr_name.fq
aligned_bam=$out_dir/$expr_name.aligned.bam
mapped_bam=$out_dir/$expr_name.aligned.mapped.bam



########### pipeline:

#### trim read. read qc
sc_trim_barcode -O $unaligned_fq -R1 $fq1 -R2 $fq2 -BS1 $index1_start -BL1 $index1_len \
    -BS2 $index2_start -BL2 $index2_len -US $umi_start -UL $umi_len -QC -N

#### alignment using subread
index_anno=test_data/ERCC92
subread-align -T 4 -i $index_anno -t 0 -r $unaligned_fq -o $aligned_bam

#### map to transcriptome
pseudo_ERCC_anno=test_data/pseudo_ERCC.gff3

# note: -BL should be inex1_len+index2_len. we dont have index1 at this experiment
# if the chr name in .bam file does not consistant with chr name in .gff3. add -C to fix the name issue
sc_exon_mapping -O $mapped_bam -I $aligned_bam -A $pseudo_ERCC_anno -BL $index2_len -UL $umi_len -S

#### demultiplexing

barcode_anno="test_data/barcode_anno.csv"

sc_demultiplex -I $mapped_bam -O $out_dir -A $barcode_anno

#### molecular counting

sc_gene_counting -O $out_dir -A $barcode_anno


