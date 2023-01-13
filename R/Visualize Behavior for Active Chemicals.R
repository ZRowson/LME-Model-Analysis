#--------------------------------------------
# Visualize chemicals with particular profiles
# Zachary Rowson
# Rowson.Zachary@epa.gov
# Created: 1/24/2022
# Last Edit: 1/24/2022
#--------------------------------------------

ac_L <- ac[1:5]
ac_D <- ac[9:13]
ac_T <- ac[14:16]
ac_trans <- ac[6:8]

hits.dt <- as.data.table(hits, keep.rownames = "cpid")
hits_long <- melt(hits.dt, id.vars = "cpid", variable.name = "acid", value.name = "hitcall")

cpid_L <- unique(hits_long[acid%in%ac_L & hitcall>0.49, cpid])
cpid_D <- unique(hits_long[acid%in%ac_D & hitcall>0.49, cpid])
cpid_T <- unique(hits_long[acid%in%ac_T & hitcall>0.49, cpid])
cpid_trans <- unique(hits_long[acid%in%ac_trans & hitcall>0.49, cpid])

cpid_LnD <- unique(c(cpid_L[which(cpid_L%in%cpid_D)]))
cpid_L.only <- cpid_L[-which(cpid_L%in%cpid_D)]
cpid_D.only <- cpid_D[-which(cpid_D%in%cpid_L)]
cpid_T.only <- cpid_T[-which(cpid_T%in%c(cpid_L,cpid_D))]
cpid_trans.only <- cpid_trans[-which(cpid_trans%in%c(cpid_L,cpid_D,cpid_T))]


library(gabi)
chm <- "Loperamide"
chm <- "5-Fluorouracil"
chm <- "5,5-diphenylhydandoin"
gabi::plot_tSeries(pmr0, chemical = chm, prsp = "SS")
chm <- cpid_L.only[1]
hits_long[cpid == chm & hitcall > 0.49, acid]
lapply(list(pmr0,pmr0_A,pmr0_J), gabi::plot_tSeries, chemical = chm)


