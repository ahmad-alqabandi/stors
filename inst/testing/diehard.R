
dh <- dieharder
summary(dh)
plot(dh)


rngs <- c("rnorm", "R_marsaglia_multic",
                    "R_super_duper", "R_mersenne_twister",
                    "R_knuth_taocp", "R_knuth_taocp2")

  rl <- lapply(rngs, function(rng) dieharder(rng, "diehard_runs", seed=12345))

 oldpar <- par(mfrow=c(2,3), mar=c(2,3,3,1))
 invisible(lapply(rl, function(res) {
  
    
   qqplot(res$data, seq(0, 1, length.out=length(res$data)),
   main=paste(res$generator, ":", round(res$p.value, digits=3)),
   ylab="", type="S")
    
    }))
 par(oldpar) # reset graph defaults



