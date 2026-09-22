data{
	int <lower=2> K;
	int n_environment;
	int xenv [K];
	int xs [K];
	int ys [K];
	int zs [K];
	matrix [n_environment,2] r;
	real <lower=0> a_beta;
	real a_b_beta;
	real b_b_beta;
	real a_overdispersion;
	real b_overdispersion;
	vector <lower=0> [K] cnv1;  // copy number for H1 (xs) per sample
	vector <lower=0> [K] cnv2;  // copy number for H2 (ys) per sample
}

parameters{
	real <lower=0> overdispersion;
	real <lower=0> b_beta;
	vector <lower=0> [K] bbeta;
	vector <lower=0> [n_environment] alpha;
}

model{
	for(k in 1:K) {
		// H1: cnv1 * (1/alpha) * r1 * beta
		xs[k] ~ neg_binomial_2(cnv1[k] * (1/alpha[xenv[k]]) * r[xenv[k],1] * bbeta[k], 1/overdispersion);
		// H2: cnv2 * alpha * r2 * beta
		ys[k] ~ neg_binomial_2(cnv2[k] * alpha[xenv[k]] * r[xenv[k],2] * bbeta[k], 1/overdispersion);
		// NonHS: reads not resolved to either haplotype, each weighted by its copy number
		zs[k] ~ neg_binomial_2((cnv1[k] * (1-r[xenv[k],1]) * (1/alpha[xenv[k]]) + cnv2[k] * (1-r[xenv[k],2]) * alpha[xenv[k]]) * bbeta[k], 1/overdispersion);
	}
	bbeta ~ gamma(a_beta,b_beta);
	b_beta ~ gamma(a_b_beta,b_b_beta);
	alpha ~ lognormal(0.0, 0.5);
	overdispersion ~ inv_gamma(a_overdispersion, b_overdispersion);
}

generated quantities{
	vector [n_environment] theta;
	for(i in 1:n_environment){
		theta[i]=1/(alpha[i]^2+1);
	}
}
