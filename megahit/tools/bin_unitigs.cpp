/*
 *  MEGAHIT
 *  Copyright (C) 2014 - 2015 The University of Hong Kong & L3 Bioinformatics Limited
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* contact: Dinghua Li <dhli@cs.hku.hk> */

#include <string>
#include <map>
#include <vector>
#include <cstdio>
#include <cassert>
#include <zlib.h>
#include <fstream>
#include <sstream>
#include "kseq.h"

using namespace std;

#ifndef KSEQ_INITED
    #define KSEQ_INITED
    KSEQ_INIT(gzFile, gzread)
#endif

char Comp(char c);
string RevComp(const string &s);
string NodeName(int i, int len, double mul, bool is_rc);

int main_bin_unitigs(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <kmer_size> <k_{kmer_size}.unitigs.fa> <genome contig probability file> <output base name> <posterior probability threshold>\n", argv[0]);
        exit(1);
    }

    unsigned k = atoi(argv[1]);
    gzFile fp = gzopen(argv[2], "r");
    assert(fp != NULL);
    kseq_t *seq = kseq_init(fp); // kseq to read files

    float pprob_threshold = atof(argv[5]);
    vector<string> ctgs;
    vector<double> muls;
    vector<string> node_names;
    vector<string> rev_node_names;
    map<string, vector<int> > start_kmer_to_id;

    while (kseq_read(seq) >= 0) {
        if (seq->seq.l < k + 1) {
            continue;
        }

        double mul;
        assert(sscanf(seq->comment.s + 7, "multi=%lf", &mul) == 1);

        muls.push_back(mul);
        ctgs.push_back(string(seq->seq.s));
    }

    for (int i = 0; i < (int)ctgs.size(); ++i) {
        start_kmer_to_id[ctgs[i].substr(0, k)].push_back(i + 1);
        start_kmer_to_id[RevComp(ctgs[i].substr(ctgs[i].length() - k))].push_back(-i - 1);
    }

    for (int i = 0; i < (int)ctgs.size(); ++i) {
        node_names.push_back(NodeName(i + 1, ctgs[i].length(), muls[i], false));
        rev_node_names.push_back(NodeName(i + 1, ctgs[i].length(), muls[i], true));
    }

    ifstream g_probs_file(argv[3]);
    string line;
    vector< vector<float> > gprobs;
    size_t genomes=0;
    while(getline(g_probs_file,line)){
        istringstream line_str(line);
        float prob;
        size_t g = 0;
        while(line_str >> prob){
            if(gprobs.size()<=g){
                gprobs.push_back(vector<float>());
                genomes++;
            }
            gprobs[g].push_back(prob);
            g++;
        }
    }

    // write out the graphs in FastG
    for(size_t g=0; g<genomes; g++){
        std::ostringstream g_fname;
        g_fname << (const char*)argv[4] << "." << g;
        ofstream g_file(g_fname.str().c_str());
        for (int i = 0; i < (int)ctgs.size(); ++i) {
            if(gprobs[g][i] < pprob_threshold) continue;
            for (int dir = 0; dir < 2; ++dir) {
                string header = dir == 0 ? node_names[i] : rev_node_names[i];
                header = ">" + header;
                string s = dir == 0 ? ctgs[i] : RevComp(ctgs[i]);
                auto mit = start_kmer_to_id.find(s.substr(s.length() - k));

                if (mit != start_kmer_to_id.end()) {
                    char sep = ':';
                    for (unsigned j = 0; j < mit->second.size(); ++j) {
                        // do not include the next node if it is below the probability threshold
                        if(gprobs[g][abs(mit->second[j])-1] < pprob_threshold) continue;
                        header += sep;
                        sep = ',';

                        if (mit->second[j] > 0) {
                            header += node_names[mit->second[j] - 1];
                        }
                        else {
                            header += rev_node_names[-mit->second[j] - 1];
                        }
                    }
                }

                header += ";";
                g_file << header << "\n" << s << "\n";
            }
        }
        g_file.close();
    }

    // write out the merged sequences
    for(size_t g=0; g<genomes; g++){
        std::ostringstream g_fname;
        g_fname << (const char*)argv[4] << ".unitigs." << g;
        ofstream g_file(g_fname.str().c_str());
        vector<bool> used(ctgs.size());
        for (int i = 0; i < (int)ctgs.size(); ++i) {
            if(gprobs[g][i] < pprob_threshold) continue;
            if(used[i]) continue;
            used[i] = true;
            string cur_seq = ctgs[i];
            ostringstream header_fwd, header_rev;
            header_fwd << i;
            for (int dir = 0; dir < 2; ++dir) {
                while(true){
                    auto mit = start_kmer_to_id.find(cur_seq.substr(cur_seq.length() - k));

                    if (mit == start_kmer_to_id.end()) break; // no outgoing edges
                    vector<int> viable;
                    for(size_t j = 0; j<mit->second.size(); j++){
                        if(gprobs[g][abs(mit->second[j])-1] >= pprob_threshold){
                            viable.push_back(mit->second[j]);
                        }
                    }

                    vector<int> rev_viable;
                    auto mit_rev = start_kmer_to_id.find(RevComp(cur_seq.substr(cur_seq.length() - k)));
                    for(size_t j = 0; j<mit_rev->second.size(); j++){
                        if(gprobs[g][abs(mit_rev->second[j])-1] >= pprob_threshold){
                            rev_viable.push_back(mit_rev->second[j]);
                        }
                    }

                    if (viable.size() != 1 || rev_viable.size() != 1) break; // need 1 outgoing edge to merge unitigs
                    
                    // extend the current sequence
                    string next_seq = ctgs[abs(viable[0])-1];
                    if(viable[0] < 0) next_seq = RevComp(next_seq);
                    cur_seq += next_seq.substr(k);
                    used[abs(viable[0])-1]=true;
                    if(dir==0) header_fwd << "," << viable[0];
                    if(dir==1) header_rev << viable[0] << ",";
                }

                cur_seq = RevComp(cur_seq);
            }
            g_file << ">" << header_rev.str() << header_fwd.str() << "\n" << cur_seq << "\n";
        }
        g_file.close();
    }
    return 0;
}
