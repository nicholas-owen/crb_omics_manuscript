#############################################################
# Description: A function to expand gene ranges
#############################################################


ExpandGeneRange <- function(x, upstream = 2000, downstream = 0)
{
    strand_is_minus = strand(x) == "-"
    on_plus = which(!strand_is_minus)
    on_minus = which(strand_is_minus)
    strand_is_minus = strand(x) == "-"
    on_plus = which(!strand_is_minus)
    on_minus = which(strand_is_minus)
    start(x)[on_plus] = start(x)[on_plus] - upstream
    start(x)[on_minus] = start(x)[on_minus] - downstream
    end(x)[on_plus] = end(x)[on_plus] + downstream
    end(x)[on_minus] = end(x)[on_minus] + upstream
    x
}
