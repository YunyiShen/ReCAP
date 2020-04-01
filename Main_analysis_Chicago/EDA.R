## plot Harvest data
Harv.data = read.csv("./_data_/Culling_1992_2008.csv",row.names = 1)

Female = colSums(Harv.data[1:8,])
Male = colSums(Harv.data[1:3+8,])

png("Harvest_data.png",width = 3.25,
    height = 2,res = 1200,units = "in",pointsize = 6)

plot(1992:2008,Female,ylab = "Harvest",xlab = "year",ylim = c(0,550))
lines(1992:2008,Female)

points(1992:2008,Male,pch = 2)
lines(1992:2008,Male,lty=2)

legend("topright", legend=c("Female","Male"),
        pch = c(0,2),lty=c(1,2), cex=1)
dev.off()









