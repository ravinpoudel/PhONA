#' Run Lasso Models
#'
#' This function takes in phyloseq object, association matix, p value matrix
#' and create a combined OTU-OTU and OTU-Phenotype network. User can select model to
#' define OTU-Phontype assocaition.
#' @param x A phyloseq object which combined OTU count, taxonomy and metadata
#' @return A matrix of the infile
#' @export

rolePlot <- function(graph.object) {
df_roles <- graph.object$role
ggplot(df_roles, aes(participation,connectivity,label=name)) + geom_point(aes(colour = factor(Treatment)))+
  scale_size(guide = 'none')+
  scale_color_manual(breaks= c("Nongraft", "Selfgraft", "RST-04-106", "Maxifort"),values = c("mediumaquamarine","sienna4","#ffae19","blueviolet"))+
  geom_text(aes(label=ifelse(role!="Peripheral",as.character(Order),'')),hjust=0.1,vjust=0.1,size=3,color="black",family="Times New Roman")+
  geom_hline(yintercept=2.5, linetype="dashed", color = "red")+
  geom_vline(xintercept = 0.62, linetype="dashed",color = "blue")+
  annotate("text", x = 0.2, y = 3, label = "Module hubs",color="gray20",size=8,family="Times New Roman")+
  annotate("text", x = 0.75, y = 3, label = " Network hubs",color="gray20",size=8,family="Times New Roman")+
  annotate("text", x = 0.2, y = -2, label = "Peripherals",color="gray20",size=8,family="Times New Roman")+
  annotate("text", x = 0.75, y = -2, label = "Connectors",color="gray20",size=8,family="Times New Roman")+
  labs(x = "Among-module connectivity", y = "Within-module degree", col="Rootstocks") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 1.0))+
  theme(axis.text=element_text(size=16,face="bold",family="Times New Roman"),axis.title=element_text(size=16,face="bold",family="Times New Roman"))+
  labs(fill = "Dose (mg)") +
  theme(legend.title=element_text(size=16,face="bold",family="Times New Roman"),
        legend.text=element_text(size=16,family="Times New Roman")) +
  guides(color = guide_legend(override.aes = list(size=5)))
}
