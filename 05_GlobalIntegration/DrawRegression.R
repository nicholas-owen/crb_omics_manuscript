
library("ggplot2")

DrawRegression <- function(x,title="")
{
    n1 <- nrow(x[[1]])
    n2 <- nrow(x[[2]])

    drawData <- data.frame(bin = c(c(1:n1),c(1:n2)), meth=c(x[[1]]$x,x[[2]]$x), tissue=c(rep("Mutation",n1),rep("WildType",n2)))
    ggplot(drawData, aes(x = bin, y = meth,  colour=tissue)) +
        geom_point(size = 2, alpha = 0.8, col=c(rep("#E69F00",n1),rep("#56B4E9",n2))) + labs(x = "average(TPM) grouped in 30 bin") +
        ylim(0, 1) +
        ylab("average(beta)") +
        geom_smooth(method = lm, se=FALSE, fullrange=TRUE) +
        scale_color_manual(values=c('#E69F00', '#56B4E9'))+
        ggtitle(title) +
        theme_minimal(base_size = 14) +
        stat_cor(label.x = 3, label.y = 0.25) +
        stat_regline_equation(label.x = 3, label.y = 0)
}

