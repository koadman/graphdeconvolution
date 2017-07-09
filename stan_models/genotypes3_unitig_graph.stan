data {
	int <lower=0> V;         // number of kmers sites
	int <lower=0> S;         // number of samples
	int <lower=0> G;         // number of genomes
	real <lower=0> depths[V*S]; // matrix of unitig depths in each sample
    int maxdepth;       // maximum depth of any genome
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
}

transformed data {
	real <lower=0> cc[V,S]; // matrix of unitig depths at each site & timepoint

    // transform the depth vector to a 2D matrix
	for(x in 1:V*S){
        int v = ((x-1)/S)+1;
        if(v <= V){
		    cc[v,((x-1)%S)+1] = depths[x];
        }
	}
}

parameters {	
    real <lower=0,upper=100> mu;        // mu parameter for lognormal strain abundance
    real <lower=0> sigma;               // sigma for lognormal strain abundance
    real <lower=0,upper=1> p_deadend;   // probability that a dead-end exists: when a node with a single outgoing vertex exists but the destination vertex does not
    real<lower=genome_min, upper=genome_max> genome_size;   // unknown average strain genome size
	vector <lower=0,upper=maxdepth>[G] abund[S];    // strain relative abundance at each time point
    vector<lower=0,upper=1>[G] geno[V];             // presence or absence of each unitig in each strain
    real<lower=0> relevance[G];     // relevance terms
}

model {
	// set the priors
    relevance ~ gamma(5,2);
	for( s in 1:S ) {
    	for( g in 1:G ) {
    		abund[s][g] ~ normal(0.0,1.0/relevance[g]); // half-normal on strain abundance terms
        }
	}
	for( v in 1:V ) {
    	for( g in 1:G ) {
		    geno[v][g] ~ beta(1.0/relevance[g],0.1); // prior concentrated on 0 & 1 (poor man's approximation to what should be binary variable)
        }
	}
    p_deadend ~ beta(0.1,0.5);

	// calculate log likelihood
	for( v in 1:V) {
		for( s in 1:S ) {
			// likelihood of Poisson distributed depths with variance-stabilizing transformation
			target += normal_log( sqrt(cc[v,s]), sqrt(abund[s]'*geno[v]), 0.25);
		}
	}

    // likelihood of nodes in each strain with 1 outgoing edge
    for( g in 1:G ){
        for( b in 1:adj1count ) {
            real s1p = geno[abs(adj1source[b])][g];
            real d1p = geno[abs(adj1dest[b])][g];
            real edge_prob = (p_deadend)*( s1p*(1-d1p) );
            edge_prob = edge_prob + (1-p_deadend)*(s1p*d1p + (1-s1p)*(1-d1p) + (1-s1p)*d1p);
            target += log(edge_prob);
        }
    }

    // likelihood of nodes in each strain with 2 outgoing edges
    for( g in 1:G ){
        for( b in 1:adj2count ) {
            real s1p = geno[abs(adj2source[b])][g];
            real d1p = geno[abs(adj2dest1[b])][g];
            real d2p = geno[abs(adj2dest2[b])][g];
            real edge_prob = (p_deadend)*( s1p*(1-d1p)*(1-d2p) );   // deadend if source exists but dests do not
            edge_prob = edge_prob + (p_deadend)*( s1p*d1p*d2p );    // hack: deadend if all three nodes exist -- this could be split out to a different term
            // not deadend if the source node does not exist, or if the source exists and only one of the dest nodes exist
            edge_prob = edge_prob + (1-p_deadend)*(s1p*d1p*(1-d2p) + s1p*(1-d1p)*d2p + (1-s1p)*(1-d1p)*(1-d2p) + (1-s1p)*d1p*(1-d2p) + (1-s1p)*(1-d1p)*d2p + (1-s1p)*d1p*d2p);
            target +=  log(edge_prob);
        }
    }

    // likelihood of genome length
    for( g in 1:G ){
        real glen = 0;
        for( v in 1:V ){
            glen = glen + geno[v][g]*lengths[v];
        }
        target +=  normal_log(glen, genome_size, genome_size/10);
    }
}

generated quantities {
    real glen[G];    
    for( g in 1:G ){
        glen[g] = 0;
        for( v in 1:V ){
            glen[g] = glen[g] + round(geno[v][g])*lengths[v];
        }
    }
}

