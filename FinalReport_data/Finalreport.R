library(dplyr)
library(reshape2)
library(ggplot2)
library(survival)
library(survminer)

######### Survival plot(PDF 5*10) ######################

survival <- read.csv("/Users/nam-yunju/Dropbox/lecture/2022-2_R_lecture/final_report/survival_plot/survival.csv")
head(survival)

s <- ggsurvplot(fit = survfit(Surv(PFS_MONTHS, PFS_STATUS)~SUBTYPE, data = survival),
                xlab="Months",
                ylab="Progression-free survival(%)",
                legend.labs=c("IDH mutation and 1p/19q codeletion", "IDH mutation and no 1p/19q codeletion", "Wild type IDH"),
                legend.title="Subtype", 
                censor=T,
                palette = c("#D5573B", "#885053", "#777DA7")) 

s$plot + theme(legend.title=element_text(size=13, face = "bold"),
               legend.text=element_text(size=13),
               legend.position = c(0.8, 0.8)) +
  scale_y_continuous(labels = scales::percent)



######### Mutationlandscape plot (PDF 10*15)######################
# only paper sample, remove no subtype info samples
mutation_data <- read.csv("/Users/nam-yunju/Desktop/mutations.txt", sep = "\t")
# contain all subtype info
subtype <- read.csv("/Users/nam-yunju/Desktop/all_subtype.csv")

head(mutation_data)
head(subtype)

colnames(subtype)[1] <- "Patient.ID"
colnames(subtype)[2] <- "Subtype"

merged <- merge(mutation_data, subtype, by="Patient.ID")
head(merged)

# filter top 19 genes 
merged_extracted <- merged %>% select(Hugo_Symbol, Variant_Classification, Patient.ID, Subtype) %>% filter(Hugo_Symbol%in%c("IDH1", "IDH2", "TP53","ATRX","CIC","NOTCH1","FUBP1","PIK3CA","NF1","PIK3R1","PTEN","EGFR","SMARCA4","ARID1A","TCF12","ZBTB20","PTPN11","PLCG1","ZCCHC12"))

# Order y-axis of plot
merged_extracted$Hugo_Symbol <-  factor(merged_extracted$Hugo_Symbol, levels = rev(c("IDH1", "IDH2", "TP53","ATRX","CIC","NOTCH1","FUBP1","PIK3CA","NF1","PIK3R1","PTEN","EGFR","SMARCA4","ARID1A","TCF12","ZBTB20","PTPN11","PLCG1","ZCCHC12")))

# Change dataframe format to order x-axis
changed <- dcast(merged_extracted, Patient.ID ~ Hugo_Symbol, value.var = "Variant_Classification")
unique(merged_extracted$Hugo_Symbol)

order_num <- order(changed[,"IDH1"], changed[,"IDH2"], changed[,"TP53"], changed[,"ATRX"], changed[,"CIC"], changed[,"NOTCH1"], changed[,"FUBP1"], changed[,"PIK3CA"], changed[,"NF1"], changed[,"PIK3R1"], changed[,"PTEN"], changed[,"EGFR"], changed[,"SMARCA4"], changed[,"ARID1A"], changed[,"TCF12"], changed[,"ZBTB20"], changed[,"PTPN11"], changed[,"PLCG1"], changed[,"ZCCHC12"])

levels <- changed$Patient.ID[order_num]
levels

merged_extracted$Patient.ID <- factor(merged_extracted$Patient.ID, levels = levels)

head(merged_extracted)

merged_extracted <- merged_extracted %>% mutate(Subtype = recode(Subtype, 
                                                                 "LGG_IDHwt" = "IDH wild type",
                                                                 "LGG_IDHmut-non-codel" = "IDH mutation, no 1p/19q codeletion",
                                                                 "LGG_IDHmut-codel" = "IDH mutation, 1p/19q codeletion"))

p <- ggplot(data=merged_extracted, aes(x=Patient.ID, y=Hugo_Symbol, fill=Variant_Classification)) + geom_tile() +
  facet_grid(cols = vars(Subtype), scales = "free", space = "free") +
  scale_fill_manual(
    values=c("Frameshift"="#647AA3","In_Frame_indel"="#CC998D", "Missense_Mutation" = "#EF798A", "Nonsense_Mutation"="#ffd6a5", "Silent"="#8e7dbe", "Splice_Site"="#76c893"),
    name="Mutation type",
    labels=c("Frameshift","In-Frame indel", "Missense", "Nonsense", "Synonymous", "Splice site")) +
  ylab("") + xlab("") +
  geom_hline(yintercept=0.5:20, colour="white") + # 가로선
  scale_y_discrete(expand=c(0,0)) + scale_x_discrete(expand=c(0,0)) + # plot내 위 아래 간격 
  theme(axis.ticks=element_line(color = NA), axis.text.x=element_blank(),
     #   axis.text=element_text(margin=unit(0.5,"cm"), colour="black", face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold", color = "black"),
        #   legend.key = element_rect(colour = NA),
        legend.text=element_text(size=12),
        legend.title=element_text(size=13, face = "bold"),
        panel.grid.major = element_blank(), # 세로선
        #   panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 17, color = "white"), # panel title
        #   strip.background = element_rect(color = "white", size = 1),
        panel.spacing = unit(0.1, "cm"), # panel사이의 간격
        legend.position = "bottom",
        plot.margin = unit(c(1, 2, 1, 2), "cm"))

g <- ggplot_gtable(ggplot_build(p))

stripr <- which(grepl("strip-t-1", g$layout$name)|grepl("strip-t-2", g$layout$name)|grepl("strip-t-3", g$layout$name))
fills <- c("#4D8C7C", "#7A8C69", "#B7D998")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)



######### arm deletion plot(PDF 10*15) ######################

arm <- read.csv("/Users/nam-yunju/Desktop/arm.csv")
unique(armlevel$Status)

arm_loss <- subset(arm, Status=="Loss")

arm_loss <- arm_loss %>% mutate(Subtype = recode(Subtype, 
                                                 "LGG_IDHwt" = "IDH wild type",
                                                 "LGG_IDHmut-non-codel" = "IDH mutation, no codeletion",
                                                 "LGG_IDHmut-codel" = "IDH mutation, codeletion"))

head(arm_loss)

a <- ggplot(data = arm_loss, aes(x=Arm)) + geom_bar(stat="count", fill = "#A38B90") +
  facet_grid(rows = vars(Subtype)) +
  ylab("Number of patients") + xlab("Arm") +
  scale_x_discrete(limits =c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q","7q","8p","8q","9p","9q","10p","10q","11p","11q","12p","12q","13q","14q","15q","16p","16q","17p","17q","18p","18q","19p","19q","20p","20q","21q","22q")) +
  theme(panel.grid.major = element_blank(),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 13),
  #      axis.text.y = element_text(size = 16, face = "bold"),
  #      axis.text.x = element_text(size = 16, face = "bold"),
        panel.spacing = unit(0.2, "cm"),
        strip.text = element_text(face = "bold", size = 15, color = "white"),
    #   panel.background = element_rect(fill = alpha("#93A39F", 0.2)),
        plot.margin = unit(c(1, 1.5, 1, 1.5), "cm")) +
  scale_y_continuous(breaks = c(0,50,100,150,200))

g_fin <- ggplot_gtable(ggplot_build(a))

stripr <- which(grepl("strip-r-1", g_fin$layout$name)|grepl("strip-r-2", g_fin$layout$name)|grepl("strip-r-3", g_fin$layout$name))
fills <- c("#4D8C7C", "#7A8C69", "#B7D998")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g_fin$grobs[[i]]$grobs[[1]]$childrenOrder))
  g_fin$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g_fin)







