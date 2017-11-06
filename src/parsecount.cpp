// parse count
#include "parsecount.h"

using std::string;

std::unordered_map<string, std::vector<std::pair<std::string,int>>> read_count(string fn, char sep)
{
    std::ifstream infile(fn);
    std::unordered_map<string, std::vector<std::pair<std::string,int>>> gene_read;
    string line;
    std::getline(infile, line); // skip header
    std::vector<string> tmp_ln;
    while(std::getline(infile, line))
    {
        std::stringstream linestream(line);
        string gene_id;
        string UMI;
        int pos;
        tmp_ln = split(line,sep);
        gene_id = tmp_ln[0];
        UMI = tmp_ln[1];
        pos = std::stoi(tmp_ln[2]);

        gene_read[gene_id].push_back(std::make_pair(UMI,pos));

    }
    infile.close();
    return gene_read;
}

int UMI_correct1(std::map<std::pair<string,int>, int>& UMI_count)
{
    bool found = false;
    int corrected_UMI = 0;
    for (auto UMI1 = UMI_count.begin(); UMI1 != UMI_count.end();) // use normal iterator in order to use `erase`
    {
        found = false;
        for (auto const& UMI2: UMI_count) // use range based
        {
            if (hamming_distance(UMI1->first.first, UMI2.first.first) == 1) // sequencing errors
            {
                if (UMI1->second == 1 || UMI1->second < UMI2.second*0.1)
                {
                    found = true;
                    // merge two UMIs
                    UMI_count[UMI2.first] += UMI_count[UMI1->first];
                    if (__DEBUG){std::cout << "merge: " <<  "<"<< UMI1->first.first << ", " << UMI1->first.second << ">" << "::" << "<"<< UMI2.first.first << ", " << UMI2.first.second << ">" << "\t" << UMI1->second << "::" << UMI2.second << std::endl;}
                    break;
                }
            } else if (UMI1->first.first == UMI2.first.first) // possiblely same sequence in different position
            {
                if ((UMI1->first.second != UMI2.first.second) && (UMI1->second == 1 || UMI1->second < UMI2.second*0.5))
                {
                    found = true;
                    // merge two UMIs
                    UMI_count[UMI2.first] += UMI_count[UMI1->first];
                    if (__DEBUG){std::cout << "merge: " <<  "<"<< UMI1->first.first << ", " << UMI1->first.second << ">" << "::" << "<"<< UMI2.first.first << ", " << UMI2.first.second << ">" << "\t" << UMI1->second << "::" << UMI2.second << std::endl;}
                    break;
                }
            }
        }
        if (found)
        {
            // delete UMI1
            corrected_UMI++;
            UMI1 = UMI_count.erase(UMI1);
        }
        else
        {
            UMI1++;
        }
    }
    return corrected_UMI;
}

int UMI_correct2(std::map<std::pair<string,int>, int>& UMI_count)
{
    bool found = false;
    int corrected_UMI = 0;
    for (auto UMI1 = UMI_count.begin(); UMI1 != UMI_count.end();) // use normal iterator in order to use `erase`
    {
        found = false;
        for (auto const& UMI2: UMI_count) // use range based
        {
            if ((UMI1->first.second - UMI2.first.second>-2 && UMI1->first.second - UMI2.first.second < 2) && (hamming_distance(UMI1->first.first, UMI2.first.first) == 1)) // sequencing errors
            {
                if (UMI1->second == 1 || UMI1->second < UMI2.second*0.5)
                {
                    found = true;
                    // merge two UMIs
                    UMI_count[UMI2.first] += UMI_count[UMI1->first];
                    if (__DEBUG){std::cout << "merge: " <<  UMI1->first.first << "::" << UMI2.first.first << "\t" << UMI1->second << "::" << UMI2.second << std::endl;}
                    break;
                }
            }else if ((UMI1->first.second - UMI2.first.second>-2 && UMI1->first.second - UMI2.first.second < 2 && UMI1->first.second != UMI2.first.second) && (UMI1->first.first == UMI2.first.first)) // diff pos
            {
                if (UMI1->second == 1 || UMI1->second < UMI2.second*0.5)
                {
                    found = true;
                    // merge two UMIs
                    UMI_count[UMI2.first] += UMI_count[UMI1->first];
                    if (__DEBUG){std::cout << "merge: " <<  UMI1->first.first << "::" << UMI2.first.first << "\t" << UMI1->second << "::" << UMI2.second << std::endl;}
                    break;
                }
            }
        }
        if (found)
        {
            // delete UMI1
            corrected_UMI++;
            UMI1 = UMI_count.erase(UMI1);
        }
        else
        {
            UMI1++;
        }
    }
    return corrected_UMI;
}


std::unordered_map<string, int> UMI_dedup(std::unordered_map<string, std::vector<std::pair<string,int>>> gene_read, std::vector<int>& UMI_dup_count, UMI_dedup_stat& s, int UMI_correct, bool read_filter)
{
    std::unordered_map<string, int> gene_counter;

    for(auto const& a_gene: gene_read)
    {
        if (read_filter && a_gene.second.size() == 1)
        {
            s.filtered_gene++;
            continue;
        }

        std::map<std::pair<std::string,int>, int> UMI_count = vector_counter(a_gene.second);
        if (UMI_correct == 1)
        {
            s.corrected_UMI += UMI_correct1(UMI_count);
        }else if (UMI_correct == 2)
        {
            s.corrected_UMI += UMI_correct2(UMI_count);
        }
        else
        {
            std::cout << "not implemented" << std::endl;
            exit(1);
        }

        for (auto const& UMI: UMI_count)
        {
            if (UMI.second > MAX_UMI_DUP)
            {
                UMI_dup_count[MAX_UMI_DUP] ++;
            }
            else
            {
                UMI_dup_count[UMI.second-1] ++;
            }

        }

        //TODO: add ATCG percentage
        gene_counter[a_gene.first] = UMI_count.size();
    }



    return gene_counter;
}

void write_mat(string fn, std::unordered_map<string, std::vector<int>> gene_cnt_matrix, std::vector<string> cellid_list)
{
    std::ofstream o_file(fn);
    //write header
    o_file << "gene_id";
    for (auto const& ce : cellid_list)
    {
        o_file << "," << ce;
    }
    o_file << std::endl;

    for (auto const& ge : gene_cnt_matrix)
    {
        o_file << ge.first;
        for (auto const& n : ge.second)
        {
            o_file << "," << n;
        }
        o_file << std::endl;
    }
    o_file.close();
}

void write_stat(string cnt_fn, string stat_fn, std::vector<int> UMI_dup_count, std::unordered_map<string, UMI_dedup_stat> UMI_dedup_stat_dict)
{
    std::ofstream cnt_file(cnt_fn);
    cnt_file << "duplication_number,count" << std::endl;
    for (int i=0; i<UMI_dup_count.size(); i++)
    {
        cnt_file << i+1 << "," << UMI_dup_count[i] << std::endl;
    }
    cnt_file.close();

    std::ofstream stat_file(stat_fn);
    stat_file << "cell_id,number_of_filtered_gene,number_of_corrected_UMI,UMI_A_percentage,UMI_T_percentage,UMI_G_percentage,UMI_C_percentage" << std::endl;
    for (auto const& n : UMI_dedup_stat_dict)
    {
        stat_file << n.first << "," << n.second.filtered_gene << "," << n.second.corrected_UMI << "," \
            << n.second.A_prop << "," << n.second.T_prop << "," << n.second.G_prop << "," << n.second.C_prop << "," << std::endl;
    }
    stat_file.close();
}

void get_counting_matrix(Barcode bar, string in_dir, int UMI_correct, bool read_filter)
{
    char sep = ',';
    std::unordered_map<string, string> cnt_files = bar.get_count_file_path(join_path(in_dir, "count"));
    std::unordered_map<string, std::vector<int>> gene_cnt_matrix; // store gene counting matrix
    std::vector<string> all_gene_list; // store all gene ids
    std::vector<int> UMI_dup_count(MAX_UMI_DUP+1, 0); // store UMI duplication statistics
    std::unordered_map<string, UMI_dedup_stat> UMI_dedup_stat_dict;
    int cell_number = bar.cellid_list.size();
    int ind = 0;
    for (auto const& ce : bar.cellid_list) // for each cell
    {
        UMI_dedup_stat_dict[ce] = {}; // init zero
        std::unordered_map<string, std::vector<std::pair<std::string,int>>> gene_read = read_count(cnt_files[ce], sep);
        std::unordered_map<string, int> gene_cnt =  UMI_dedup(gene_read, UMI_dup_count, UMI_dedup_stat_dict[ce], UMI_correct, read_filter);

        for (auto const& ge : gene_cnt) // for each gene
        {
            if (gene_cnt_matrix.find(ge.first) == gene_cnt_matrix.end())
            {
                auto & vec = gene_cnt_matrix[ge.first];
                vec.resize(cell_number, 0); // init with all zeros
                vec[ind] = ge.second;
            }
            else
            {
                gene_cnt_matrix[ge.first][ind] = ge.second;
            }
        }
        ind++;
    }

    // write to file
    write_mat(join_path(in_dir, "gene_count.csv"), gene_cnt_matrix, bar.cellid_list);
    string stat_dir = join_path(in_dir, "stat");
    write_stat(join_path(stat_dir, "UMI_duplication_count.csv"), join_path(stat_dir, "UMI_dedup_stat.csv"), UMI_dup_count, UMI_dedup_stat_dict);

}

