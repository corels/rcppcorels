
#include "queue.h"

#define STRICT_R_HEADERS
#include <Rcpp.h>

#define BUFSZ 512

/*
 * Logs statistics about the execution of the algorithm and dumps it to a file.
 * To turn off, pass verbosity <= 1
 */
NullLogger* logger;

//' R Interface to 'Certifiably Optimal RulE ListS (Corels)'
//'
//' Nothing here yet.
//' @title Corels interace
//' @param rules_file Character variable with rules file name
//' @param labels_file Character variable with labels file name
//' @param log_dir Character variable with logfile directory name
//' @param meta_file Optional character variable with labels file name
//' @param run_bfs Boolean toggle
//' @param calculate_size Boolean toggle
//' @param run_curiosity Boolean toggle
//' @param curiosity_policy Integer value
//' @param latex_out Boolean toggle
//' @param map_type Integer value
//' @param verbosity Integer value
//' @param max_num_nodes Integer value
//' @param regularization Double value
//' @param logging_frequency Integer value
//' @param ablation Double value
//' @return A constant bool for now
//' @examples
//' library(RcppCorels)
//'
//' logdir <- tempdir()
//' rules_file <- system.file("sample_data", "compas_train.out", package="RcppCorels")
//' labels_file <- system.file("sample_data", "compas_train.label", package="RcppCorels")
//' meta_file <- system.file("sample_data", "compas_train.minor", package="RcppCorels")
//'
//' stopifnot(file.exists(rules_file),
//'           file.exists(labels_file),
//'           file.exists(meta_file),
//'           dir.exists(logdir))
//'
//' corels(rules_file, labels_file, logdir, meta_file,
//'        verbosity = 100,
//'        regularization = 0.015,
//'        curiosity_policy = 2,   # by lower bound
//'        map_type = 1) 	   # permutation map
//' cat("See ", logdir, " for result file.")
// [[Rcpp::export]]
bool corels(std::string rules_file,
            std::string labels_file,
            std::string log_dir,
            std::string meta_file = "",
            bool run_bfs = false,
            bool calculate_size = false,
            bool run_curiosity = false,
            int curiosity_policy = 0,
            bool latex_out = false,
            int map_type = 0,
            int verbosity = 0,
            int max_num_nodes = 100000,
            double regularization = 0.01,
            int logging_frequency = 1000,
            int ablation = 0) {

    //bool run_bfs = false;
    //bool run_curiosity = false;
    //int curiosity_policy = 0;
    //bool latex_out = false;
    bool use_prefix_perm_map = (map_type==1);
    bool use_captured_sym_map = (map_type==2);
    //int verbosity = 0;
    //int map_type = 0;
    //int max_num_nodes = 100000;
    double c = /*0.01*/ regularization;
    //char ch;
    //bool error = false;
    //char error_txt[BUFSZ];
    int freq = /*1000*/ logging_frequency;
    //int ablation = 0;
    //bool calculate_size = false;


    std::map<int, std::string> curiosity_map;
    curiosity_map[1] = "curiosity";
    curiosity_map[2] = "curious_lb";
    curiosity_map[3] = "curious_obj";
    curiosity_map[4] = "dfs";

    int nrules, nsamples, nlabels, nsamples_chk;
    rule_t *rules, *labels;
    rules_init(rules_file.c_str(), &nrules, &nsamples, &rules, 1);
    rules_init(labels_file.c_str(), &nlabels, &nsamples_chk, &labels, 0);

    int nmeta, nsamples_check;
    // Equivalent points information is precomputed, read in from file, and stored in meta
    rule_t *meta;
    if (meta_file != "")
        rules_init(meta_file.c_str(), &nmeta, &nsamples_check, &meta, 0);
    else
        meta = NULL;

    if (verbosity >= 10)
        print_machine_info();
    char froot[BUFSZ];
    char log_fname[BUFSZ+4];
    char opt_fname[BUFSZ+8];
    const char* pch = strrchr(rules_file.c_str(), '/');
    snprintf(froot, BUFSZ, "%s/for-%s-%s%s-%s-%s-removed=%s-max_num_nodes=%d-c=%.7f-v=%d-f=%d",
             log_dir.c_str(),
             pch ? pch + 1 : "",
             run_bfs ? "bfs" : "",
             run_curiosity ? curiosity_map[curiosity_policy].c_str() : "",
             use_prefix_perm_map ? "with_prefix_perm_map" :
                (use_captured_sym_map ? "with_captured_symmetry_map" : "no_pmap"),
             meta ? "minor" : "no_minor",
             ablation ? ((ablation == 1) ? "support" : "lookahead") : "none",
             max_num_nodes, c, verbosity, freq);
    snprintf(log_fname, BUFSZ+4, "%s.txt", froot);
    snprintf(opt_fname, BUFSZ+8, "%s-opt.txt", froot);

    if (verbosity >= 1000) {
        Rprintf("\n%d rules %d samples\n\n", nrules, nsamples);
        rule_print_all(rules, nrules, nsamples);

        Rprintf("\nLabels (%d) for %d samples\n\n", nlabels, nsamples);
        rule_print_all(labels, nlabels, nsamples);
    }

    if (verbosity > 1)
        logger = new Logger(c, nrules, verbosity, log_fname, freq);
    else
        logger = new NullLogger();
    double init = timestamp();
    char run_type[BUFSZ];
    Queue* q;
    strcpy(run_type, "LEARNING RULE LIST via ");
    char const *type = "node";
    if (curiosity_policy == 1) {
        strcat(run_type, "CURIOUS");
        q = new Queue(curious_cmp, run_type);
        type = "curious";
    } else if (curiosity_policy == 2) {
        strcat(run_type, "LOWER BOUND");
        q = new Queue(lb_cmp, run_type);
    } else if (curiosity_policy == 3) {
        strcat(run_type, "OBJECTIVE");
        q = new Queue(objective_cmp, run_type);
    } else if (curiosity_policy == 4) {
        strcat(run_type, "DFS");
        q = new Queue(dfs_cmp, run_type);
    } else {
        strcat(run_type, "BFS");
        q = new Queue(base_cmp, run_type);
    }

    PermutationMap* p;
    if (use_prefix_perm_map) {
        strcat(run_type, " Prefix Map\n");
        PrefixPermutationMap* prefix_pmap = new PrefixPermutationMap;
        p = (PermutationMap*) prefix_pmap;
    } else if (use_captured_sym_map) {
        strcat(run_type, " Captured Symmetry Map\n");
        CapturedPermutationMap* cap_pmap = new CapturedPermutationMap;
        p = (PermutationMap*) cap_pmap;
    } else {
        strcat(run_type, " No Permutation Map\n");
        NullPermutationMap* null_pmap = new NullPermutationMap;
        p = (PermutationMap*) null_pmap;
    }

    CacheTree* tree = new CacheTree(nsamples, nrules, c, rules, labels, meta, ablation, calculate_size, type);
    Rprintf("%s", run_type);
    // runs our algorithm
    bbound(tree, max_num_nodes, q, p);

    Rprintf("final num_nodes: %zu\n", tree->num_nodes());
    Rprintf("final num_evaluated: %zu\n", tree->num_evaluated());
    Rprintf("final min_objective: %1.5f\n", tree->min_objective());
    const tracking_vector<unsigned short, DataStruct::Tree>& r_list = tree->opt_rulelist();
    Rprintf("final accuracy: %1.5f\n",
       1 - tree->min_objective() + c*r_list.size());
    print_final_rulelist(r_list, tree->opt_predictions(),
                     latex_out, rules, labels, opt_fname);

    Rprintf("final total time: %f\n", time_diff(init));
    logger->dumpState();
    logger->closeFile();
    if (meta) {
        Rprintf("\ndelete identical points indicator");
        rules_free(meta, nmeta, 0);
    }
    Rprintf("\ndelete rules\n");
    rules_free(rules, nrules, 1);
    Rprintf("delete labels\n");
    rules_free(labels, nlabels, 0);
    Rprintf("tree destructors\n");


    return true;                  // more to fill in, naturally
}
