data {
	int <lower=0> V;         // number of kmers sites
	int <lower=0> S;         // number of samples
	int <lower=0> G;         // number of genomes
	real <lower=0> depths[V*S]; // matrix of unitig depths in each sample
    int maxdepth;   // maximum depth of any genome
    int adj1count;
    int adj1source[adj1count];
    int adj1dest[adj1count];
    int adj2count;
    int adj2source[adj2count];
    int adj2dest1[adj2count];
    int adj2dest2[adj2count];
    real lengths[V];
    real expected_length;
}

transformed data {
	real <lower=0> cc[V,S]; // matrix of unitig depths at each site & timepoint

	for(x in 1:V*S){
        int v = ((x-1)/S)+1;
        if(v <= V){
		    cc[v,((x-1)%S)+1] = depths[x];
        }
	}
}

parameters {	
    real <lower=0,upper=100> mu;
    real <lower=0> sigma;
    real <lower=0,upper=1> p_deadend; // probability that a dead-end exists: when a node with a single outgoing vertex exists but the destination vertex does not
    real<lower=0, upper=10000> genome_size;
	vector <lower=0,upper=maxdepth>[G] abund[S]; // strain relative abundance at each time point
    vector<lower=0,upper=1>[G] geno[V];
}

model {
	// set the priors
	for( s in 1:S ) {
		abund[s] ~ lognormal(mu,sigma); // lognormally distributed strain abundances
	}
	for( v in 1:V ) {
    	for( g in 1:G ) {
		    geno[v][g] ~ beta(0.1,0.1); // prior concentrated on 0 & 1 (poor man's approximation to what should be binary variable)
        }
	}
    genome_size ~ normal(1400,50);
    p_deadend ~ beta(0.1,0.5);

	// calculate log likelihood
	for( v in 1:V) {
		for( s in 1:S ) {
			// Poisson distributed depths with variance-stabilizing transformation
			target += normal_log( sqrt(cc[v,s]), sqrt(abund[s]'*geno[v]), 0.25);
		}
	}

    // likelihood of graph structure
    for( g in 1:G ){
        for( b in 1:adj1count ) {
            real s1p = geno[abs(adj1source[b])][g];
            real d1p = geno[abs(adj1dest[b])][g];
            real edge_prob = (p_deadend)*( s1p*(1-d1p) );
            edge_prob = edge_prob + (1-p_deadend)*(s1p*d1p + (1-s1p)*(1-d1p) + (1-s1p)*d1p);
            target += log(edge_prob);
        }
    }

    for( g in 1:G ){
        for( b in 1:adj2count ) {
            real s1p = geno[abs(adj2source[b])][g];
            real d1p = geno[abs(adj2dest1[b])][g];
            real d2p = geno[abs(adj2dest2[b])][g];
            real edge_prob = (p_deadend)*( s1p*(1-d1p)*(1-d2p) );
            edge_prob = edge_prob + (p_deadend)*( s1p*d1p*d2p );
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
        target +=  normal_log(glen, expected_length, expected_length/30);
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

