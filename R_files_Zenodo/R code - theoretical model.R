


f.soln <- function(t.prior1, t.rr.pred, t.rr.none, t.s1, t.u1, t.u2, t.k1, t.k2){
	
	p.b = t.prior1	
		
	if(!isTRUE(all.equal(0, t.u1))){
		mu.result = uniroot(f.mu, c(-10,10), uu = t.u1)$root
		mu.a1 = -mu.result
		mu.b1 = mu.result
	} else {mu.a1 = -999; mu.b1 = 999}
	
		
	if(!isTRUE(all.equal(0, t.u2))){
		mu.result2 = uniroot(f.mu, c(-10,10), uu = t.u2)$root
		mu.a2 = -mu.result2
		mu.b2 = mu.result2
	}else {mu.a2 = -999; mu.b2 = 999}

	
	if(t.rr.pred == 0 && t.rr.none == 0){
		beta = (1-p.b)/p.b
		s1c = (2*(log(p.b/(1-p.b))) - mu.b1^2 + mu.a1^2) / (2*(mu.a1 - mu.b1)) 
	}	else {
		beta = (1-p.b)*t.rr.none/(p.b*t.rr.pred)
		s1c = (2*(log(p.b*t.rr.pred/((1-p.b)*t.rr.none))) - mu.b1^2 + mu.a1^2) / (2*(mu.a1 - mu.b1))
		}

	
	
	p.hit = 1 - pnorm(s1c, mu.b1, sigma1)
	p.fa  = 1 - pnorm(s1c, mu.a1, sigma1)


	f.p.1c = function(x){ 
	# function for solving for cut-off probability p.1c in the absence of stimuli
		x*(1-p.b)/(p.b*(1-x)) - beta
	}

	p.1c = ifelse(beta == Inf, 1,
			ifelse(beta == 0, 0, uniroot(f.p.1c, c(0, 1))$root)) # solve for p.1c, the cutoff prob in the absence of stimuli and account for the limits of beta
	
	
	vi.1.actual = ifelse(p.b >= p.1c,
		p.b*p.hit*t.rr.pred - (1-p.b)*p.fa*t.rr.none + p.b*t.rr.pred*(1-beta) - t.k1,
		p.b*p.hit*t.rr.pred - (1-p.b)*p.fa*t.rr.none - t.k1
	)
	vi.1 = vi.1.actual
	int1 = vi.1 >= 0 # is the first stimulus used?


	########
	# update prior and new cutoff stimulus if first stimulus is used
	if(int1 == T ){
		p.s1 = p.b*dnorm(t.s1, mu.b1, sigma1) + (1-p.b)*dnorm(t.s1, mu.a1, sigma1)
		p.b2 = p.b*dnorm(t.s1, mu.b1, sigma1) / p.s1

	}
		
	#######
	# don't update prior if first stimulus is ignored
	if(int1 == F){
		p.b2 = p.b # prior isn't updated because 1st stimulus was ignored
	}
	

	########
	# 2nd stimulus

	if(t.rr.pred == 0 && t.rr.none == 0){
		beta2 = (1-p.b2)/p.b2
		s2c = (2*( log( p.b2/(1-p.b2) )) - mu.b2^2 + mu.a2^2) / (2*(mu.a2 - mu.b2)) 
	}	else {
		beta2 = (1-p.b2)*t.rr.none/(p.b2*t.rr.pred)
		s2c = (2*( log( p.b2*t.rr.pred/((1-p.b2)*t.rr.none) )) - mu.b2^2 + mu.a2^2) / (2*(mu.a2 - mu.b2))
		}
	
	p.hit2 = 1 - pnorm(s2c, mu.b2, sigma2)
	p.fa2 = 1 - pnorm(s2c, mu.a2, sigma2)	
	

	f.p.2c = function(x){
		x*(1-p.b2)/(p.b2*(1-x)) - beta2
	}
	
	p.2c = ifelse(beta2 == Inf, 1, 
			ifelse(isTRUE(all.equal(0, beta2)), 0, uniroot(f.p.2c, c(0, 1))$root))  
	
	vi.2.actual = ifelse(p.b2 >= p.2c,
			p.b2*p.hit2*t.rr.pred - (1-p.b2)*p.fa2*t.rr.none + p.b2*t.rr.pred*(1-beta2) - t.k2,
			p.b2*p.hit2*t.rr.pred - (1-p.b2)*p.fa2*t.rr.none - t.k2
		)
		
	vi.2 = ifelse(isTRUE(all.equal(0, vi.2.actual)), 0, vi.2.actual)
	
	int2 = vi.2 >= 0 

	
	####### which stimulus/stimuli was/were used? (0 = neither, 1 = 1st only, 2 = 2nd only, 3 = both)
	st.use = ifelse (int1 == T && int2 == T, 3, # integration
					ifelse (int1 == T && int2 == F, 1, # only 1st stimulus used
						ifelse (int1 == F && int2 == T, 2, # only 2nd stimulus used
							0))) # neither stimulus was used
		
	list(int1 = int1, vi.1 = vi.1, int2 = int2, vi.2 = vi.2, st.use = st.use)

}