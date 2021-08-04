library(trialr)
# ps <- list()
# ps[[1]] <- list(p.true=c(0.25, 0.40, 0.5, 0.55, 0.6),
#                 pE.true=c(0.15, 0.25, 0.3, 0.45, 0.55))
# ps[[2]] <- list(p.true=c(0.25, 0.40, 0.5, 0.55, 0.6),
#                 pE.true=c(0.15, 0.25, 0.5, 0.5, 0.5))
# ps[[3]] <- list(p.true=c(0.25, 0.40, 0.5, 0.55, 0.6),
#                 pE.true=c(0.15, 0.25, 0.5, 0.3, 0.1))
# ps[[4]] <- list(p.true=c(0.1, 0.15, 0.25, 0.30, 0.45),
#                 pE.true=c(0.05, 0.10, 0.17, 0.25, 0.35))
# ps[[5]] <- list(p.true=c(0.22, 0.45, 0.55, 0.65, 0.70),
#                 pE.true=c(0.03, 0.10, 0.20, 0.35, 0.40))
#  
# # 2
# idx <- 2
# p.true <- ps[[idx]]$p.true
# pE.true <- ps[[idx]]$pE.true

res.efftox.fn <- function(res, ndose=5){
    nsim  <- length(res$recommended_dose)
    Sel <- sapply(1:ndose, function(i){sum(res$recommended_dose==i, na.rm=TRUE)})/nsim
    Sel <- round(Sel*100, 1)
    Allo <- sapply(1:ndose, function(i){sum(unlist(res$doses_given)==i)})/nsim
    ntoxs <- sum(unlist(res$toxicities))/nsim
    neffs <- sum(unlist(res$efficacies))/nsim
    nsubs <- sum(Allo)
    errStop <- 100 - sum(Sel)
                
    sum.res.efftox <- list(Selection=Sel,
                           Allocation=Allo,
                           errStop=errStop,
                           tol.effs=neffs,
                           tol.toxs=ntoxs,
                           to.Subjs=nsubs)
    return(sum.res.efftox)
}

eff.p <- efftox_solve_p(eff0 = 0.5, tox1 = 0.65, eff_star = 0.7, tox_star = 0.25)
eff.dat <- list(num_doses = 5,
            real_doses = c(1, 2, 4, 6.6, 10),
            efficacy_hurdle = phiE,
            toxicity_hurdle = phi,
            p_e = 0.1,
            p_t = 0.1,
            p = eff.p,
            eff0 = 0.5,
            tox1 = 0.65,
            eff_star = 0.7,
            tox_star = 0.25,
            alpha_mean = -7.9593, alpha_sd = 3.5487,
            beta_mean = 1.5482, beta_sd = 3.5018,
            gamma_mean = 0.7367, gamma_sd = 2.5423,
            zeta_mean = 3.4181, zeta_sd = 2.4406,
            eta_mean = 0, eta_sd = 0.2,
            psi_mean = 0, psi_sd = 1,
            doses = c(),
            tox   = c(),
            eff   = c(),
            num_patients = 0)

res.efftox <- efftox_simulate(eff.dat, 
                              num_sims=nsimu, 
                              first_dose=1, 
                              true_eff=pE.true,
                              true_tox=p.true,
                              cohort_sizes=rep(cohortsize, ncohort), 
                              chains=4, cores=4)

sum.res.efftox <- res.efftox.fn(res.efftox)
sum.res.efftox$p.true <- p.true
sum.res.efftox$pE.true <- pE.true
