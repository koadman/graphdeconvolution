data {
	int <lower=0,upper=1>	using_reads;	// switch between read count and kmer count models
	int <lower=0> V;         // number of unitigs
	int <lower=0> S;         // number of samples
	int <lower=0> G;         // number of genomes
	int <lower=0> readcounts[V,S]; // matrix of mapped read counts for unitigs in each sample
	real <lower=0> kmercov[V,S]; // matrix of unitig coverage in each sample
  int adj1count;      // number of nodes with outdegree 1
  int adj1source[adj1count];  // source node with outdegree 1
  int adj1dest[adj1count];    // dest node w/outdegree 1
  int adj2count;      // number of nodes with outdegree 2
  int adj2source[adj2count];  // source node w/outdegree 2
  int adj2dest1[adj2count];   // first dest node w/outdegree 2
  int adj2dest2[adj2count];   // second dest node w/outdegree 2
  real lengths[V];        // length in nucleotides of each node
  real<lower=0> genome_min;   // smallest genome size
  real<lower=0> genome_max;   // largest genome size
	real <lower=0,upper=1> p_deadend;   // probability that a dead-end exists: when a node with a single outgoing vertex exists but the destination vertex does not
	real <lower=0,upper=1> p_fork;   // probability that the graph forks due to a repeat sequence in a strain
	int read_len;

	// data for link model
	int L; // number of links among sites (edges in graph)
	int links[L,2]; // IDs of unitigs that have been co-observed in read pairs
	int linkcounts[L,S]; // counts of unitig co-observation, at each time point

	int coretigcount;
	int coretigs[coretigcount];

}

parameters {
    real<lower=genome_min, upper=genome_max> genome_size;   // unknown average strain genome size
    vector<lower=0,upper=1>[G] geno[V];             // presence or absence of each unitig in each strain
    real<lower=0,upper=1> relevance[G];     // relevance terms
		vector <lower=0>[G] abund[S];    // strain relative abundance at each time point
}

model {
	for( v in 1:V ) {
    for( g in 1:G ) {
			geno[v][g] ~ beta(0.01,0.01); // prior concentrated on 0 & 1 (poor man's approximation to what should be binary variable)
    }
	}

	for(ctig in 1:coretigcount){
		for( g in 1:G ) {
			geno[coretigs[ctig]][g] ~ beta(1,0.01); // core contigs are in every strain
    }
	}

	// calculate log likelihood
	for( v in 1:V) {
		for( s in 1:S ) {
			real cov = abund[s]'*geno[v]; // ' transposes the vector
			target += using_reads*lengths[v]*poisson_lpmf( readcounts[v,s] | (cov*(lengths[v]+read_len)) / read_len);
			// likelihood of Poisson distributed depths with variance-stabilizing transformation
			target += (1-using_reads)*lengths[v]*normal_lpdf( sqrt(kmercov[v,s]) | sqrt(cov), 0.25);
		}
	}
	// likelihood that a node is due to sequencing error
//	for( v in 1:V) {
//		real g_sum;
//		g_sum = lengths[v]*sum(depths[v])*log(p_seqerror);
//	}


    // likelihood of nodes in each strain with 1 outgoing edge
    for( g in 1:G ){
        for( b in 1:adj1count ) {
            real s1p = geno[abs(adj1source[b])][g];
            real d1p = geno[abs(adj1dest[b])][g];
            real edge_prob = (p_deadend)*( s1p*(1-d1p) );
            edge_prob = edge_prob + (1-p_deadend)*(s1p*d1p + (1-s1p)*(1-d1p) + (1-s1p)*d1p);
//            target += log(edge_prob);
        }
    }

    // likelihood of nodes in each strain with 2 outgoing edges
    for( g in 1:G ){
        for( b in 1:adj2count ) {
            real s1p = geno[abs(adj2source[b])][g];
            real d1p = geno[abs(adj2dest1[b])][g];
            real d2p = geno[abs(adj2dest2[b])][g];
            real edge_prob = (p_deadend)*( s1p*(1-d1p)*(1-d2p) );   // deadend if source exists but dests do not
            edge_prob = edge_prob + (p_fork)*( s1p*d1p*d2p );    // fork if all three nodes exist
						edge_prob = edge_prob + (p_deadend)*( (1-s1p)*d1p*(1-d2p) + (1-s1p)*(1-d1p)*d2p );
						edge_prob = edge_prob + (p_deadend*p_deadend)*( (1-s1p)*d1p*d2p );
            // no deadends if the none of the nodes exist, or if the source exists and only one of the dest nodes exist
            edge_prob = edge_prob + (1-p_deadend)*(s1p*d1p*(1-d2p) + s1p*(1-d1p)*d2p + (1-s1p)*(1-d1p)*(1-d2p));
//            target +=  log(edge_prob);
        }
    }

    // likelihood of genome length
    for( g in 1:G ){
        real glen = 0;
        for( v in 1:V ){
            glen = glen + geno[v][g]*lengths[v];
        }
        target +=  normal_lpdf(glen | genome_size, genome_size/10);
    }

		// likelihood of link configurations
		for( i in 1:L ){
			for( s in 1:S ){
				vector[G] p_haves;
				for( g in 1:G ){
					// for each link, joint probability that:
					// 1. strains lack one or both in the pair and all links are erroneous
					p_haves[g] = geno[links[i,1]][g]*geno[links[i,2]][g];
				}
				// number of expected links between the node pair, based upon
				// the presence or absence of nodes in strains and the strain abundances
				// FIXME this needs to incorporate the read fragment size distribution
//				target += poisson_lpmf( linkcounts[i,s] | abund[s]'*p_haves);
			}
    }
}
