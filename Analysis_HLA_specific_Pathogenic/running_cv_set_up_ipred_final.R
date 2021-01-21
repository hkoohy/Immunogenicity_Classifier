
library(data.table)
library(dplyr)
library(stringdist)
library(reshape2)
library(parallel)
library(EMCluster)
library(ggplot2)
library(pROC)
select <- dplyr::select
CORES <- 4
res_em=NULL
# Substitute chowell training data for new training data. For simplicity, keep the name as 'chowell'
chowell <- fread("classifier/GS_iPred_train.txt")
# kidera below can also be calculated by using compute kidera function below for each amino acid.
kidera <- fread("classifier/kidera.txt") %>% melt


compute_kidera <- function(data) {
  if (is.data.frame(data)) {
    if ("antigen.epitope" %in% colnames(data)) {
      data <- data$antigen.epitope
    } else {
      stop("The 'antigen.epitope' column is missing in data frame")
    }
  } else {
    if (!is.character(data)) {
      stop("Unknown input type, should be either character vector or a data frame with 'antigen.epitope' column")
    }
  }
  res <- data %>%
    strsplit("") %>%
    mclapply(function(x) data.table(antigen.epitope = paste0(x, collapse = ""),
                                    aa = x) %>%
               merge(kidera, allow.cartesian = T) %>%
               group_by(antigen.epitope, variable) %>%
               summarise(value = sum(value)),
             mc.cores = CORES) %>%
    rbindlist %>%
    dcast(antigen.epitope ~ variable, fun.aggregate = sum)
}

# test
compute_kidera(c("A"))

chowell.k <- merge(chowell, compute_kidera(chowell))

mat.chowell.k <- chowell.k %>% 
  select(-antigen.epitope, -hla, -immunogenicity) %>%
  as.matrix
pc <- prcomp(mat.chowell.k, 
             scale = T, rank = 2)
chowell.k.pc <- chowell.k
chowell.k.pc$pcX <- pc$x[,1]
chowell.k.pc$pcY <- pc$x[,2]



pcaplot=ggplot(chowell.k.pc, aes(x = pcX, y = pcY)) +
  stat_density_2d(data = chowell.k.pc %>% select(pcX, pcY), geom = "raster",
                  aes(fill = ..density..), contour = F) +
  geom_density2d(aes(color = immunogenicity)) +
  scale_color_brewer("Immunogenic", palette = "Set1", labels = c("No", "Yes")) +
  scale_fill_gradient(low = "white", high="black") +
  scale_x_continuous(expand=c(0,0), limits = c(-4,4))+
  scale_y_continuous(expand=c(0,0), limits = c(-4,4))+
  theme_bw() + 
  theme(aspect = 1,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +labs(fill="Density") + xlab("PC-2") + ylab("PC-1") + font("xy.text",size=16,color="black")+
  font("xlab",size=16,color="black")+ font("ylab",size=16,color="black")  + font("legend.text",size=14) + font("legend.title",size=14)
print(pcaplot)

mat.chowell.k

res_em <- init.EM(mat.chowell.k, 
                  nclass = 2,
                  lab = ifelse(chowell.k$immunogenicity == "Positive", 1, 2))
summary(res_em)
res_chowell <- e.step(mat.chowell.k, res_em) %>% 
  as.data.table %>%
  cbind(chowell.k)

densityplot=res_chowell %>%
  ggplot(aes(x = Gamma.V1, fill = immunogenicity)) +
  geom_density(alpha = 0.9) +
  scale_fill_brewer("Immunogenic", palette = "Set1", labels = c("No", "Yes")) +
  scale_x_continuous("P(Immunogenic)") +
  scale_y_continuous("Density") +
  theme_bw() +
  theme(aspect = 1, 
        legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +ggtitle(paste0("AUC: ",roc(immunogenicity ~ Gamma.V1, res_chowell, ci = F, plot = F)$auc))
print(densityplot)


roc(immunogenicity ~ Gamma.V1, res_chowell, ci = F, plot = F)

predict_imm <- function(data) {
  if (is.data.frame(data)) {
    data.k <- merge(data, compute_kidera(data))
  } else {
    data.k <- compute_kidera(data)
  }
  mat.k <- data.k %>% select(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10) %>% as.matrix
  data.k$imm.prob <- e.step(mat.k, res_em) %>% as.data.table %>% .$Gamma.V1
  data.k %>% select(-f1, -f2, -f3, -f4, -f5, -f6, -f7, -f8, -f9, -f10)
}

predict_imm(c("ELALGIGILV", "LAPGATNEK", "GLCTVAML"))

