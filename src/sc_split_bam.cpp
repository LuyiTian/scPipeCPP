#include "parsebam.h"
#include "cellbarcode.h"
#include "timing.h"

using std::string;

int main(int argc, char* argv[]) {
    if (argc < 6) 
    { // Insufficient parameters
        std::cout << "Usage: sc_split_bam \n" <<\
            "\t-I <infile> the input bam file (required)\n"<<\
            "\t-O <outdir> the output dir (required)\n"<<\
            "\t-A <barcode_annotation_file> annotate barcode and cell_id, \n for Drop-seq data you should run `./sc_detect_bc` to get the annotation file (required)\n"<<\
            "\t-M <max_mismatch> maximum mismatch allowed when matching barcode (default: 1)\n"<<\
            "\t-T <max_tmp_reads> maximum number of temporal reads to keep in memory (default: 1000000)\n"<<\
            "\t-BC <cellular_tag> two characters cell barcode tag used in bam file (default: `XC`)\n";
        exit(0);
    } 
    else 
    {
        string bam_fn, out_dir, anno_fn;
        string bc = "YC";
        string mb = "YM";
        string gb = "GE";
        string am = "YE";
        string mt = "chrM";
        int max_mismatch = 1;
        int max_tmp_reads = 1000000;
        bool has_UMI = true;
        for (int i = 1; i < argc; i++) 
        {
            string arg = argv[i];

            if (i + 1 != argc) 
            {
                if (arg == "-I") 
                {
                    bam_fn = argv[i + 1];
                } 
                else if (arg == "-A") 
                {
                    anno_fn = argv[i + 1];
                }
                else if (arg == "-O") 
                {
                    out_dir = argv[i + 1];
                } 
                else if (arg == "-M")
                {
                    max_mismatch = argv[i + 1][0]-'0';
                }
                else if (arg == "-T")
                {
                    max_tmp_reads = std::atoi(argv[i + 1]);
                }
                else if (arg == "-BC")
                {
                    bc = argv[i + 1];
                }
            }
        }
        std::cout << "######### demultiplexing:" << std::endl;
        std::cout << "parameters:" << std::endl;
        std::cout << "\t" << "bam file: " << bam_fn << std::endl;
        std::cout << "\t" << "out dir: " << out_dir << std::endl;
        std::cout << "\t" << "annotation file: " << anno_fn << std::endl;
        std::cout << "\t" << "max mismatch: " << max_mismatch << std::endl;
        std::cout << "\t" << "max temporal reads: " << max_tmp_reads << std::endl;
        std::cout << "\t" << "cell tag: " << bc << std::endl;

        Timer timer;
        timer.start();
        timer.mode("minutes");

        Barcode bar;
        bar.read_anno(anno_fn);
        Bamdemultiplex bam_de = Bamdemultiplex(out_dir, bar, bc, mb, gb, am, mt);
        bam_de.split_bam(bam_fn, max_mismatch, max_tmp_reads);

        std::cout << "Time elapsed: " << timer.end() << std::endl;
        
        return 0;
    }
}