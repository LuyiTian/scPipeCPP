// sc_exon_mapping.cpp

#include "transcriptmapping.h"
#include "timing.h"

int main(int argc, char* argv[]) {
    if (argc < 4)
    { // Insufficient parameters
        std::cout << "Usage: sc_exon_mapping \n" <<\
            "\t-O <outbam> the output bam file (required)\n"<<\
            "\t-I <inbam> aligned bam file (required)\n"<<\
            "\t-A <exon_annotation> gff3 or bed annotation files separated with comma.(required)\n"<<\
            "\t-AM <map_status_tag> mapping status tag (default: YE)\n"<<\
            "\t-GE <gene_id_tag> gene id tag (default: GE)\n"<<\
            "\t-BC <barcode_tag> cell barcode tag (default: YC)\n"<<\
            "\t-MB <molecular_barcode> molecular barcode tag (default: YM)\n"<<\
            "\t-BL <barcode_length> total barcode length  (required)\n"<<\
            "\t-UL <UMI_length> UMI length  (required)\n"<<\
            "\t-S <strand_specific> perform strand specific mapping (default: false)\n"<<\
            "\t-C <chrname_fix> add `chr` to chromosome names, `MT` to `chrM` (default: false)\n";
        exit(0);
    }
    else
    {
        std::string anno_files, fn, fn_out, am, ge, bc, mb;
        bc = "YC";
        mb = "YM";
        ge = "GE";
        am = "YE";
        int bc_len = -1;
        int UMI_len = -1;
        bool m_strand = false;
        bool fix_chr = false;
        for (int i = 1; i < argc; i++)
        {
            std::string arg = argv[i];

            if (i + 1 != argc)
            {
                if (arg == "-O")
                {
                    fn_out = argv[i + 1];
                }
                else if (arg == "-I")
                {
                    fn = argv[i + 1];
                }
                else if (arg == "-A")
                {
                    anno_files = argv[i + 1];
                }
                else if (arg == "-AM")
                {
                    am = argv[i + 1];
                }
                else if (arg == "-GE")
                {
                    ge = argv[i + 1];
                }
                else if (arg == "-BC")
                {
                    bc = argv[i + 1];
                }
                else if (arg == "-MB")
                {
                    mb = argv[i + 1];
                }
                else if (arg == "-BL")
                {
                    bc_len = std::atoi(argv[i + 1]);
                }
                else if (arg == "-UL")
                {
                    UMI_len = std::atoi(argv[i + 1]);
                }
            }
            if (arg == "-S")
            {
                m_strand = true;
            }
            else if (arg == "-C")
            {
                fix_chr = true;
            }
        }
        

        std::vector<std::string> token = split(anno_files, ',');

        std::cout << "######### transcriptome mapping:" << std::endl;
        std::cout << "parameters:" << std::endl;
        std::cout << "\t" << "output bam file: " << fn_out << std::endl;
        std::cout << "\t" << "aligned input bam file file: " << fn << std::endl;
        for( const auto& n : token)
        {
            std::cout << "\t" << "annotation file: " << n << std::endl;
        }
        std::cout << "\t" << "map status tag: " << am << std::endl;
        std::cout << "\t" << "geneid tag: " << ge << std::endl;
        std::cout << "\t" << "cell barcode tag: " << bc << std::endl;
        std::cout << "\t" << "molecular barcode tag: " << mb << std::endl;
        std::cout << "\t" << "barcode length: " << bc_len << std::endl;
        std::cout << "\t" << "UMI_length: " << UMI_len << std::endl;
        std::cout << "\t" << "strand_specific: " << std::boolalpha << m_strand << std::endl;
        std::cout << "\t" << "add `chr` to chromosome names: " << std::boolalpha << fix_chr << std::endl;

        if (UMI_len<0 || bc_len<0)
        {
            std::cout << "ERROR: UMI length and barcode length must be specified." << std::endl;
            exit(1);
        }
        Mapping a = Mapping();
        for (const auto& n : token)
        {
            a.add_annotation(n, fix_chr);
        }
        Timer timer;
        timer.start();
        timer.mode("seconds");
        a.parse_align(fn, fn_out, m_strand, am, ge, bc, mb, bc_len, UMI_len);
        
        std::cout << "Time elapsed: " << timer.end() << std::endl;
        
        return 0;
    }
}