//trim_barcode
#include "trimbarcode.h"

void kseq_t_to_bam_t(kseq_t *seq, bam1_t *b, int trim_n)
{
    int seq_l = seq->seq.l - trim_n; // seq length after trim the barcode
    b->l_data = seq->name.l + 1 + (int)(1.5 * seq_l + (seq_l % 2 != 0)); // +1 includes the tailing '\0'
    if (b->m_data < b->l_data) 
    {
        b->m_data = b->l_data;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
    }
    b->core.tid = -1;
    b->core.pos = -1;
    b->core.mtid = -1;
    b->core.mpos = -1;
    b->core.flag = BAM_FUNMAP;
    b->core.l_qname = seq->name.l + 1; // +1 includes the tailing '\0'
    b->core.l_qseq = seq_l;
    b->core.n_cigar = 0; // we have no cigar sequence
    memcpy(b->data, seq->name.s, sizeof(char) * seq->name.l); // first set qname
    b->data[seq->name.l] = '\0';
    uint8_t *s = bam_get_seq(b);
    int i = 0;
    for (i = 0; i < b->core.l_qseq;++i) // set sequence
    {
        bam1_seq_seti(s, i, seq_nt16_table[seq->seq.s[i + trim_n]]);
    }

    s = bam_get_qual(b);

    for (i = 0; i < b->core.l_qseq;++i) // set quality
    {
        s[i] = seq->qual.s[i + trim_n]-33;
    }
}

static bool check_qual(char *qual_s, int trim_n, int thr, int below_thr)
{
    int not_pass = 0;
    for (int i = 0; i < trim_n; i++)
    {
        if ((int)qual_s[i] <= thr){
            not_pass++;
        }
    }
    return not_pass>below_thr?false:true;
}

static bool N_check(char *seq, int trim_n){
    bool pass = true;
    char *ptr = strchr(seq, 'N');
    if (ptr)
    {
        int index = ptr - seq;
        if (index <= trim_n)
        {
            pass = false;
        }
    }
    return pass;
}

void paired_fastq_to_bam(char *fq1_fn, char *fq2_fn, char *bam_out, const read_s read_structure, const filter_s filter_settings)
{
    // Filter tallies
    int passed_reads = 0;
    int removed_have_N = 0;
    int removed_low_qual = 0;

    int l1 = 0;
    int l2 = 0;

    gzFile fq1 = gzopen(fq1_fn, "r"); // input fastq
    if (!fq1){fprintf(stderr, "cant open file: %s\n", fq1_fn); exit(EXIT_FAILURE);}
    gzFile fq2 = gzopen(fq2_fn, "r");
    if (!fq2){fprintf(stderr, "cant open file: %s\n", fq2_fn); exit(EXIT_FAILURE);}

    samFile *fp = sam_open(bam_out,"wb"); // output file

    // write header
    bam_hdr_t *hdr = bam_hdr_init();
    hdr->l_text = strlen(empty_header);
    hdr->text = strdup(empty_header);
    hdr->n_targets = 0;
    sam_hdr_write(fp, hdr);

    // extract settings
    int id1_st = read_structure.id1_st;
    int id1_len = read_structure.id1_len;
    int id2_st = read_structure.id2_st;
    int id2_len = read_structure.id2_len;
    int umi_st = read_structure.umi_st;
    int umi_len = read_structure.umi_len;

    int bc1_end, bc2_end; // get total length of index + UMI for read1 and read2
    int state; // 0 for two index with umi, 1 for two index without umi, 2 for one index with umi, 3 for one index without umi

    const int TWO_INDEX_WITH_UMI = 0;
    const int TWO_INDEX_NO_UMI = 1;
    const int ONE_INDEX_WITH_UMI = 2;
    const int ONE_INDEX_NO_UMI = 3;

    if (id1_st >= 0) // if we have plate index
    {
        state = TWO_INDEX_WITH_UMI;
        bc1_end = id1_st + id1_len;
    }
    else // if no plate information, use id1_len to trim the read 1
    {
        state = ONE_INDEX_WITH_UMI;
        bc1_end = id1_len;
        uint8_t *idx = new uint8_t[id2_len + 1];
    }

    // set barcode end index
    if (umi_st >= 0)
    {
        if (id2_st + id2_len > umi_st + umi_len)
        {
            bc2_end = id2_st + id2_len;
        }
        else
        {
            bc2_end = umi_st + umi_len;
        }
    }
    else
    {
        state++; // no umi
        bc2_end = id2_st + id2_len;
    }

    // set offset for fastq header
    int name_offset;
    if (state == TWO_INDEX_WITH_UMI)
    {
        name_offset = id1_len + id2_len + umi_len + 2;
    }
    else if (state == TWO_INDEX_NO_UMI)
    {
        name_offset = id1_len + id2_len + 2;
    }
    else if (state == ONE_INDEX_WITH_UMI)
    {
        name_offset = id2_len + umi_len + 2;
    }
    else if (state == ONE_INDEX_NO_UMI)
    {
        name_offset = id2_len + 2;
    }

    // initialise fastq readers
    kseq_t *seq1;
    seq1 = kseq_init(fq1);
    kseq_t *seq2;
    seq2 = kseq_init(fq2);

    // main loop, iterate through each fastq record
    // assume there are the name number of reads in read1 and read2 files, not checked.
    while (((l1 = kseq_read(seq1)) >= 0) && ((l2 = kseq_read(seq2)) >= 0))
    {  
        // qual check before we do anything
        if (filter_settings.if_check_qual)
        {
            if(!(check_qual(seq1->seq.s, bc1_end, filter_settings.min_qual, filter_settings.num_below_min) && check_qual(seq2->seq.s, bc2_end, filter_settings.min_qual, filter_settings.num_below_min)))
            {
                removed_low_qual++;
                continue;
            }
        }
        if (filter_settings.if_remove_N)
        {
            if(!(N_check(seq1->seq.s, bc1_end) && N_check(seq2->seq.s, bc2_end)))
            {
                removed_have_N++;
                continue;
            }
        }
        
        // begin processing valid read
        passed_reads++;

        bam1_t *b = bam_init1();
        
        // move original read name
        int new_name_length = name_offset + seq1->name.l;
        char* old_name_adr = seq1->name.s;  
        char* new_name_adr = old_name_adr + name_offset;
        int n_char_copied = seq1->name.l * sizeof(char);
        seq1->name.s = (char*)realloc(old_name_adr, new_name_length);
        memcpy(new_name_adr, old_name_adr, n_char_copied);

        if (state == TWO_INDEX_WITH_UMI)
        {
            memcpy(seq1->name.s, seq1->seq.s + id1_st, id1_len*sizeof(char)); // copy index one
            memcpy(seq1->name.s + id1_len, seq2->seq.s + id2_st, id2_len*sizeof(char)); // copy index two
            seq1->name.s[id1_len + id2_len] = '_'; // add separator
            memcpy(seq1->name.s + id1_len + id2_len + 1, seq2->seq.s + umi_st, umi_len*sizeof(char)); // copy umi

        }
        else if (state == TWO_INDEX_NO_UMI)
        {
            memcpy(seq1->name.s, seq1->seq.s + id1_st, id1_len*sizeof(char)); // copy index one
            memcpy(seq1->name.s + id1_len, seq2->seq.s + id2_st, id2_len*sizeof(char)); // copy index two
            seq1->name.s[id1_len + id2_len] = '_'; // add separator
        }
        else if (state == ONE_INDEX_WITH_UMI)
        {
            memcpy(seq1->name.s, seq2->seq.s + id2_st, id2_len*sizeof(char)); // copy index two
            seq1->name.s[id2_len] = '_'; // add separator
            memcpy(seq1->name.s + id2_len + 1, seq2->seq.s + umi_st, umi_len*sizeof(char)); // copy umi
        }
        else if (state == ONE_INDEX_NO_UMI)
        {
            memcpy(seq1->name.s, seq2->seq.s + id2_st, id2_len*sizeof(char)); // copy index two
            seq1->name.s[id2_len] = '_'; // add separator
        }
        seq1->name.s[name_offset-1] = '#';
        seq1->name.l = name_offset + seq1->name.l;

        kseq_t_to_bam_t(seq1, b, bc1_end);
        int ret = sam_write1(fp, hdr, b);
        if (ret < 0)
        {
            std::cout << "fail to write the bam file: " << seq1->name.s << std::endl;
            std::cout << "return code: " << ret << std::endl;
            exit(EXIT_FAILURE);
        }
        bam_destroy1(b);
    }

    kseq_destroy(seq1); kseq_destroy(seq2); // free seq 
    gzclose(fq1); gzclose(fq2); // close fastq file
    sam_close(fp); // close bam file
    std::cout << "pass QC: " << passed_reads << std::endl;
    std::cout << "removed_have_N: " << removed_have_N << std::endl;
    std::cout << "removed_low_qual: " << removed_low_qual << std::endl;
}