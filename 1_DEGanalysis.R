

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


DATA <- read.csv("Gene_expression_DATA.csv") # arranged mice_data, no other difference
rownames(DATA)<-DATA[,1]
DATA<-DATA[,-1]
DATA<-as.matrix(DATA)


seqrebTAM<-seq(from=1,to=7,by=2)
reb_TAM<-DATA[,seqrebTAM]
seqrebTC<-seq(from=2,to=8,by=2)
reb_TC<-DATA[,seqrebTC]
seqepTAM<-seq(from=9,to=20,by=2)
ep_TAM<-DATA[,seqepTAM]
seqepTC<-seq(from=10,to=20,by=2)
ep_TC<-DATA[,seqepTC]
seqvehTAM<-seq(from=21,to=30,by=2)
veh_TAM<-DATA[,seqvehTAM]
seqvehTC<-seq(from=22,to=31,by=2)
veh_TC<-DATA[,seqvehTC]

mean_vTAM<-apply(veh_TAM,1,mean)
mean_vTC<-apply(veh_TC,1,mean)

log_rTAM<-log2(reb_TAM)/log2(mean_vTAM)
log_eTAM<-log2(ep_TAM)/log2(mean_vTAM)
log_rTC<-log2(reb_TC)/log2(mean_vTC)
log_eTC<-log2(ep_TC)/log2(mean_vTC)


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


##### t test in TAM-associated genes
t_reb_TAM<-log_rTAM
t_ep_TAM<-log_eTAM

P_value=matrix(NA,1,nrow(t_reb_TAM))
i_TAM<-c() 
for (i in 1:nrow(t_reb_TAM))#start t_TAM
{if(all(t_reb_TAM[i,]==t_reb_TAM[i,1])==FALSE&&
    all(t_ep_TAM[i,]==t_ep_TAM[i,1])==FALSE&&
    is.finite(t_reb_TAM[i,1])&&is.finite(t_reb_TAM[i,2])&&is.finite(t_reb_TAM[i,3])&&is.finite(t_reb_TAM[i,4])&&
    is.finite(t_ep_TAM[i,1])&&is.finite(t_ep_TAM[i,2])&&is.finite(t_ep_TAM[i,3])&&is.finite(t_ep_TAM[i,4])&&is.finite(t_ep_TAM[i,5])&&is.finite(t_ep_TAM[i,6]))
{P_value[i]<-t.test(t_reb_TAM[i,],t_ep_TAM[i,])$p.value
if(P_value[i]<0.05)
{i_TAM<-c(i_TAM,i)
}
}
}

P_adjust=p.adjust(P_value[i_TAM],method="fdr",n=length(P_value[i_TAM]))
i_TAM=i_TAM[P_adjust<0.05]
length(i_TAM) # 2869  0.01 1176

sel_t_rTAM<-t_reb_TAM[i_TAM,]
sel_t_eTAM<-t_ep_TAM[i_TAM,]


##### fold change 
mean_rTAM<-apply(sel_t_rTAM,1,mean)
mean_eTAM<-apply(sel_t_eTAM,1,mean)

fc_TAM<-abs(log2(mean_rTAM/mean_eTAM))


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


##### t test in TC-associated genes
t_reb_TC<-log_rTC
t_ep_TC<-log_eTC

P_value=matrix(0,1,nrow(t_reb_TC))
i_TC<-c() 
#pv<-c()
for (i in 1:nrow(t_reb_TC))#start t_TC
{if(all(t_reb_TC[i,]==t_reb_TC[i,1])==FALSE&&
    all(t_ep_TC[i,]==t_ep_TC[i,1])==FALSE&&
    is.finite(t_reb_TC[i,1])&&is.finite(t_reb_TC[i,2])&&is.finite(t_reb_TC[i,3])&&is.finite(t_reb_TC[i,4])&&
    is.finite(t_ep_TC[i,1])&&is.finite(t_ep_TC[i,2])&&is.finite(t_ep_TC[i,3])&&is.finite(t_ep_TC[i,4])&&is.finite(t_ep_TC[i,5])&&is.finite(t_ep_TC[i,6]))
{P_value[i]<-t.test(t_reb_TC[i,],t_ep_TC[i,])$p.value
if(P_value[i]<0.05)
{i_TC<-c(i_TC,i)
}
}
}

P_adjust=p.adjust(P_value[i_TC],method="fdr",n=length(P_value[i_TC]))
i_TC=i_TC[P_adjust<0.05]
length(i_TC) #6512   0.01 3398

sel_t_rTC<-t_reb_TC[i_TC,]
sel_t_eTC<-t_ep_TC[i_TC,]


##### fold change 
mean_rTC<-apply(sel_t_rTC,1,mean)
mean_eTC<-apply(sel_t_eTC,1,mean)

fc_TC<-abs(log2(mean_rTC/mean_eTC))


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


namesTAM=names(fc_TAM[fc_TAM>log2(1.5)]) # 310

saveRDS(namesTAM, "mice_namesTAM.rds")

humanTAM=toupper(namesTAM)

#####
Marker_genes=humanTAM
Marker_genes=intersect(Marker_genes,rownames(Gene_CGGALGG))
Marker_genes=intersect(Marker_genes,rownames(Gene_LGG))

mapTAM=humanTAM[-c(match(Marker_genes,humanTAM))]

orTAM=match(mapTAM,humanTAM)

map=cbind(orTAM,mapTAM)
#####

humanTAM[20]=c('MTFR2')#
humanTAM[24]=c('CD24')
humanTAM[44]=c('NLRP1')
humanTAM[62]=c('SEPSECS')
humanTAM[68]=c('C14orf143')#
humanTAM[82]=c('ARHGEF28')#
humanTAM[88]=c('SAMD4A')
humanTAM[97]=c('C8orf55')#
humanTAM[122]=c('HLA-DQA1')
humanTAM[124]=c('C4A')
humanTAM[142]=c('MYRF')#
humanTAM[188]=c('FAM212B')#
humanTAM[197]=c('ZNF462')
humanTAM[258]=c('TENM4')#
humanTAM[263]=c('SEPP1')

write.csv(humanTAM,file="humanTAM.csv")


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


namesTC=names(fc_TC[fc_TC>log2(1.5)])    # 920

saveRDS(namesTC, "mice_namesTC.rds")

humanTC=toupper(namesTC)

#####
Marker_genes=humanTC
Marker_genes=intersect(Marker_genes,rownames(Gene_CGGAGBM))
Marker_genes=intersect(Marker_genes,rownames(Gene_GBM))

mapTC=humanTC[-c(match(Marker_genes,humanTC))]

orTC=match(mapTC,humanTC)

map=cbind(orTC,mapTC)
#####

humanTC[12]=c('C2orf40')
humanTC[45]=c('NOS1AP')
humanTC[47]=c('C1orf192')
humanTC[69]=c('KIAA1244')
humanTC[98]=c('NT5DC3')
humnaTC[134]=c('ZNF692')
humanTC[172]=c('C17orf58')
humanTC[182]=c('PRCD')
humanTC[190]=c('C2orf44')
humanTC[207]=c('SPTB')
humanTC[220]=c('AKR1E2')
humanTC[222]=c('ZNF184')
humanTC[247]=c('SERF1A')
humanTC[257]=c('C14orf37')
humanTC[259]=c('SKA3')
humanTC[265]=c('PCP2')
humanTC[292]=c('ZNF623')
humanTC[301]=c('KIAA1644')
humanTC[306]=c('ZNF641')
humanTC[313]=c('ZNF174')
humanTC[330]=c('KIAA1524')
humantC[337]=c('NOTCH2')
humanTC[385]=c('PCDHB19P')
humanTC[386]=c('PCDHB14')
humanTC[408]=c('C18orf54')
humanTC[412]=c('C18orf25')
humanTC[418]=c('C11orf95')
humanTC[481]=c('C20orf112')
humanTC[499]=c('CA13')
humanTC[526]=c('CA14')
humanTC[527]=c('HIST2H3D')
humanTC[529]=c('ZNF697')
humanTC[544]=c('UGT8')
humanTC[552]=c('C1orf173')
humanTC[564]=c('LPPR1')
humanTC[565]=c('ZNF462')
humanTC[612]=c('C1orf159')
humanTC[628]=c('KIAA1239')
humanTC[632]=c('KIAA1211')
humanTC[643]=c('ZNF605')
humanTC[680]=c('KIAA1549')
humanTC[686]=c('DNAH6')
humanTC[698]=c('ZNF25')
humanTC[710]=c('C12orf69')
humanTC[720]=c('ZNF583')
humanTC[726]=c('ZNF235')
humanTC[727]=c('ZNF428')
humanTC[729]=c('MIA')
humanTC[730]=c('SPTBN4')
humanTC[733]=c('C19orf33')
humanTC[777]=c('ZNF689')
humanTC[780]=c('GLI4')
humanTC[797]=c('SGK223')
humanTC[817]=c('ZNF276')
humanTC[822]=c('KIAA1377')
humanTC[823]=c('SP5')
humanTC[825]=c('ZNF558')
humanTC[828]=c('LPPR2')
humanTC[845]=c('C11orf87')
humanTC[859]=c('CA12')
humanTC[895]=c('KIF4A')
humanTC[901]=c('ZNF711')
humanTC[909]=c('SPCS3')

write.csv(humanTC,file="humanTC.csv")


