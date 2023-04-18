### WriteOutPval.R
### function that rounds p-values for nice table generation

WriteOutPval <- function(pv.v,round.min=3,round.max=300){

  round.v <- round.min:round.max
  th.v <- 10^(-round.v);

  outpv.v <- vector(length=length(pv.v));
  done.idx <- which(pv.v >= th.v[1]);

  outpv.v[done.idx] <- round(pv.v[done.idx],round.min);

  todo.idx <- setdiff(1:length(pv.v),done.idx);


  for(i in todo.idx){
    if(length(which(th.v <= pv.v[i]))>0){
     outpv.v[i] <- round(pv.v[i],round.v[min(which(th.v <= pv.v[i]))]);
    }
    else{
     outpv.v[i] <- 0;
    }
  }

  return(outpv.v);
}
