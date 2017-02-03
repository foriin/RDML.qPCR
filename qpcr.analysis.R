# install.packages("RDML")
# install.packages("reshape2")
library(RDML)
library(reshape2)
library(dplyr)
library(ggplot2)

deriv <- function(x, y){
  dy.dx <- rep(NA, length(y) - 1)
  for (i in 2:length(y)){
    dy.dx[i - 1] <- ((y[i] - y[i - 1])/(x[i] - x[i - 1]) + (y[i + 1] - y[i])/(x[i + 1] - x[i]))/2
  }
  return(dy.dx)
}

# rRNA28s <- RDML$new("~/IMG/qPCR/PIWI_RIP/2017-01-31 17-24-27_piwi-rip_28S.rdml")
rdmlob <- RDML$new("~/IMG/qPCR/PIWI_RIP/2017-02-01 12-23-56_Piwi-RIP_Copia.rdml")

snames <- data.frame(fdata.name = c("IP Control 1",
                            "IP Control 2",
                            "IP ActD 1",
                            "IP ActD 2",
                            "Mock Control 1",
                            "Mock Control 2",
                            "Mock ActD 1",
                            "Mock ActD 2",
                            "Input Control 1",
                            "Input Control 2",
                            "Input ActD 1",
                            "Input ActD 2"),
                     sample = as.character(1:12))

# tab <- rRNA28s$AsTable(
#   # Custom name pattern 'position~sample~sample.type~target~dye'
#   name.pattern = paste(
#     react$position,
#     react$sample$id,
#     private$.sample[[react$sample$id]]$type$value,
#     data$tar$id,
#     target[[data$tar$id]]$dyeId$id,
#     sep = "~"),
#   # Custom column 'quantity' - starting quantity of added sample 
#   quantity = {
#     value  <- sample[[react$sample$id]]$quantity$value
#     if (is.null(value) || is.na(value)) NULL
#     else value
#   }
# )

atab <- rdmlob$AsTable()
atab <- merge(snames, atab[, -1],  all = T)
tab <- atab[, c(2, 1, 3:ncol(atab))]
tab$fdata.name <- as.character(tab$fdata.name)
tab$fdata.name[21:22] <- c("RT-")

fdata <- rdmlob$GetFData(tab,
                          # long table format for usage with ggplot2
                          long.table = TRUE)

mdata <- rdmlob$GetFData(tab, dp.type = "mdp")

mdata <- as.data.frame(cbind(mdata$tmp[2:nrow(mdata)], apply(mdata[, 2:15], 2, function(x){
  -deriv(mdata$tmp, x)
})))
names(mdata)[1] <- "tmp"

mp <- melt(mdata, id.vars = "tmp")


png("amplification.28S.png", width = 800, height = 600)
  ggplot(fdata, aes(cyc, fluor)) +
    geom_line(aes(group = fdata.name,
                  color = fdata.name))+
    labs(title = "Amplification",
         x = "Cycle",
         y = "RFU",
         color = "id")
dev.off()


png("melting.28s.png", width = 800, height = 600)
  ggplot(mp, aes(tmp, value, col=variable))+
    geom_line()+
    labs(title = "Melting curve (derivative)",
         x = "Temperature",
         y = "-dFU/dT",
         col="id")
dev.off()
