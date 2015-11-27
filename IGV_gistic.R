#setwd('/home/jun/Programs/GISTIC2.0/melanoma/score.gistic/')

filenames <- list.files(path = ".", pattern = '^score_*', all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
for (i in filenames){
  
  a <- read.delim(i)
  
  index = data.frame(Type = rep('Amp',2), Chromosome = rep(23,2), Start = c(10000001, 30000001),
                     End = c(30000001, 40000001), X.log10.q.value. = rep(0,2), G.score = c(2,4),
                     average.amplitude = rep(0,2), frequency = rep(0,2))
  
  b <- rbind.data.frame(a, index)
  c <- b[,c(1:4,6,5,7,8)]
  c$G.score <- round(2^round(c$G.score, digits = 1)-1,2)
  
  write.table(c, paste('igv', i, sep = ''), quote = FALSE, sep = "\t",
              row.names = FALSE, append = FALSE,
              col.names = c('Type', 'Chromosome', 'Start',
                            'End', 'G-score',	'-log10(q-value)', 'Average-amplitude',	'frequency')
              )
  
}


braf <- read.delim('scores_braf.gistic')

index = data.frame(Type = rep('Amp',2), Chromosome = rep(23,2), Start = c(1000000, 3000000),
                   End = c(2000000, 4000000), X.log10.q.value. = rep(0,2), G.score = c(1,2),
                   average.amplitude = rep(0,2), frequency = rep(0,2))

braf <- rbind.data.frame(braf, index)
summary(index)
summary(braf)
colnames(braf)
braf <- braf[,c(1:4,6,5,7,8)]

head(braf)

tail(braf)
braf$G.score <- round(2^round(braf$G.score, digits = 1)-1,2)

write.table(braf, 'braf.gistic', quote = FALSE, sep = "\t ",
            row.names = FALSE,
            col.names = TRUE)