//test_simulator.cpp
#include <string>
#include <iostream>
#include "../fq_simulator.h"


const std::string pseudo_ERCCfn = "test_data/pseudo_ERCC.gff3";
const std::string barfn = "test_data/barcode_anno.csv";
const std::string ERCCfafn = "test_data/ERCC92.fa";

// copy & paste !
void print_result(bool is_good, int &passed_test, int &failed_test){
    if(is_good)
    {
        std::cout << " passed" << std::endl;
        passed_test++;
    }
    else
    {
        std::cout << " failed" << std::endl;
        failed_test++;
    }
}


bool test_FaReader()
{
    std::cout << "\ttest class FaReader ...";
    int cnt = 0;
    bool cnt_res, seq_res;
    std::string resultseq170 = "AAGTGAGGCA";

    FaReader fareader = FaReader(ERCCfafn);
    while (fareader.readone())
    {
        cnt ++;
        if (fareader.fa.name == "ERCC-00170")
        {
            seq_res = (resultseq170 == fareader.fa.seq.substr(45, 10));
        }
    }
    cnt_res = (cnt == 92);
    return cnt_res & seq_res;
}


bool test_gamma_count_matrix()
{
    std::cout << "\ttest CountSimulator::gamma_count_matrix() ...";
    CountSimulator cnt_sim(666);  //use fixed seed to reproduce data.
    std::vector<std::string> gene_v = {"aa", "bb", "cc"};
    std::vector<std::string> cell_v = {"aa", "bb", "cc", "dd"};
    cnt_sim.init_mat(gene_v, cell_v);
    cnt_sim.gamma_count_matrix(0.3, 100);
    return cnt_sim.cnt_mat(1, 1) == 42;
}


bool test_get_transcript_seq()
{
    std::cout << "\ttest Celseq2Simulator::get_transcript_seq() ...";
    Celseq2Simulator base_sim = Celseq2Simulator(pseudo_ERCCfn, barfn);
    FaReader fareader = FaReader(ERCCfafn);
    std::string testseq33, testseq116;  // test exon split in ERCC-00033 and ERCC-00116
    bool res33, res116;
    std::string resultseq33 = "TTTTTGTTTGAAAATGAGAA";
    std::string resultseq116 = "GTGACGAACGTTGTCGCTCA";
    while (fareader.readone())
    {
        for (auto const ge : base_sim.Anno.gene_dict[fareader.fa.name])
        {
            if (fareader.fa.name == "ERCC-00033")
            {
                testseq33 = base_sim.get_transcript_seq(ge, fareader.fa).substr(1012, 20);
                res33 = (testseq33 == resultseq33);
            }
            else if (fareader.fa.name == "ERCC-00116")
            {
                testseq116 = base_sim.get_transcript_seq(ge, fareader.fa).substr(981, 20);
                res116 = (testseq116 == resultseq116);
            }
        }
    }
    return res33 & res116;
}


bool test_gen_random_seq()
{
    std::cout << "\ttest Celseq2Simulator::gen_random_seq() ...";
    Celseq2Simulator base_sim = Celseq2Simulator(pseudo_ERCCfn, barfn, 666);  //use fixed seed to reproduce data.
    return base_sim.gen_random_seq(6) == "AGTGGG";  // is it platform independent? 
}


void generate_celseq2_data()
{
    std::string R1fn = "celseq2/simu_celseq2_R1.fq";
    std::string R2fn = "celseq2/simu_celseq2_R2.fq";
    Celseq2Simulator sim(pseudo_ERCCfn, barfn, 666);
    std::vector<double> param = {1, 100};
    sim.gen_gene_expression("gamma_random", param);
    sim.makefq(R1fn, R2fn, ERCCfafn);
}




int main(int argc, char* argv[])
{
    int passed_test = 0, failed_test = 0;
    bool test_result = false;

    // test Celseq2Simulator::get_transcript_seq()
    test_result = test_get_transcript_seq();
    print_result(test_result, passed_test, failed_test);

    // test FaReader()
    test_result = test_FaReader();
    print_result(test_result, passed_test, failed_test);

    // test Celseq2Simulator::gen_random_seq()
    test_result = test_gen_random_seq();
    print_result(test_result, passed_test, failed_test);

    // test CountSimulator::gamma_count_matrix()
    test_result = test_gamma_count_matrix();
    print_result(test_result, passed_test, failed_test);

    generate_celseq2_data();

    // print all test result
    std::cout << "passed_test: " << passed_test << "; " 
    << "failed_test: " << failed_test << std::endl;

    
    return 0;
}