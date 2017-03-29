
#include <Rcpp.h>
#include "trimbarcode.h"

read_s get_read_structure(Rcpp::NumericVector bs1,
                          Rcpp::NumericVector bl1,
                          Rcpp::NumericVector bs2,
                          Rcpp::NumericVector bl2,
                          Rcpp::NumericVector us,
                          Rcpp::NumericVector ul){
  read_s s = {};
  s.id1_st = Rcpp::as<int>(bs1);  // id1 start
  s.id1_len = Rcpp::as<int>(bl1);    // id1 length
  s.id2_st = Rcpp::as<int>(bs2);     // id2 start
  s.id2_len = Rcpp::as<int>(bl2);    // id2 length
  s.umi_st = Rcpp::as<int>(us);     // umi start
  s.umi_len = Rcpp::as<int>(ul);    // umi length
  return s;
}

filter_s get_filter_structure(Rcpp::NumericVector rmlow,
                              Rcpp::NumericVector rmN,
                              Rcpp::NumericVector minq,
                              Rcpp::NumericVector numbq){
  int c_rmlow = Rcpp::as<int>(rmlow);
  int c_rmN = Rcpp::as<int>(rmN);
  filter_s f = {};
  f.if_check_qual = c_rmlow==1?true:false;
  f.if_remove_N = c_rmN==1?true:false;

  f.min_qual = Rcpp::as<int>(minq);
  f.num_below_min = Rcpp::as<int>(numbq);
  return f;
}


// [[Rcpp::export]]

void rcpp_sc_trim_barcode_paired(Rcpp::CharacterVector outfq,
                                 Rcpp::CharacterVector r1,
                                 Rcpp::CharacterVector r2,
                                 Rcpp::NumericVector bs1,
                                 Rcpp::NumericVector bl1,
                                 Rcpp::NumericVector bs2,
                                 Rcpp::NumericVector bl2,
                                 Rcpp::NumericVector us,
                                 Rcpp::NumericVector ul,
                                 Rcpp::NumericVector rmlow,
                                 Rcpp::NumericVector rmN,
                                 Rcpp::NumericVector minq,
                                 Rcpp::NumericVector numbq) {

  std::string c_outfq = Rcpp::as<std::string>(outfq);
  std::string c_r1 = Rcpp::as<std::string>(r1);
  std::string c_r2 = Rcpp::as<std::string>(r2);
  read_s s = get_read_structure(bs1, bl1, bs2, bl2, us, ul);
  filter_s fl = get_filter_structure(rmlow, rmN, minq, numbq);

  paired_fastq_to_fastq((char *)c_r1.c_str(), (char *)c_r2.c_str(), (char *)c_outfq.c_str(), s, fl);
}
