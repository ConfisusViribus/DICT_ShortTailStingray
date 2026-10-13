#### POPAN Open population model on STSR by year ####

library(RMark)
library(tidyverse)
#detach("package:marked", unload=TRUE)

STSR_data <- STSR_ID_data %>%
  filter(Date >= "2014-01-01") %>%
  select(STSR_ID, Date)  %>%
  dplyr::rename(Individual_ID = STSR_ID) %>%
  dplyr::rename(date = Date) %>%
  drop_na(Individual_ID)

data = STSR_data
# Create temporal occasions based on months
data <- data %>%
  mutate(
    year = year(date),
    month = month(date),
    occasion = (year - min(year))
  )

# Create capture histories
ch_data <- data %>%
  group_by(Individual_ID) %>%
  summarize(
    first_capture = min(occasion),
    captures = list(occasion)
  )

# Convert to binary capture history
max_occasion <- max(data$occasion)
capture_histories <- matrix(0, nrow = nrow(ch_data), ncol = max_occasion)

for(i in 1:nrow(ch_data)) {
  capture_histories[i, unlist(ch_data$captures[i])] <- 1
}

# Create capture history strings
ch_strings <- apply(capture_histories, 1, paste, collapse = "")

ch_df <- data.frame(
  ch = ch_strings
)

hw.proc = RMark::process.data(ch_df, model="POPAN")
hw.ddl = RMark::make.design.data(hw.proc)

# survival process
Phi.ct = list(formula=~1) # constant
Phi.time = list(formula=~time) # year effect
# detection process
p.ct = list(formula=~1) # constant
p.time = list(formula=~time) # year effect
# entry probablity process
pent.ct = list(formula=~1)
pent.time = list(formula=~time)


# constant survival, constant recapture
Model.1 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.ct,p=p.ct, pent=pent.ct, N=N.ct),output = T,delete=F)
# constant survival, time-dependent recapture
Model.2 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.ct,p=p.time, pent=pent.ct, N=N.ct),output = T,delete=F)
# time-dependent survival, constant recapture
Model.3 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.time,p=p.ct, pent=pent.ct, N=N.ct),output = T,delete=F)
# time-dependent survival, time-dependent recapture
Model.4 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.time,p=p.time, pent=pent.ct, N=N.ct),output = T,delete=F)

# constant survival, constant recapture
Model.5 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.ct,p=p.ct, pent=pent.time, N=N.ct),output = T,delete=F)
# constant survival, time-dependent recapture
Model.6 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.ct,p=p.time, pent=pent.time, N=N.ct),output = T,delete=F)
# time-dependent survival, constant recapture
Model.7 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.time,p=p.ct, pent=pent.time, N=N.ct),output = T,delete=F)
# time-dependent survival, time-dependent recapture
Model.8 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.time,p=p.time, pent=pent.time, N=N.ct),output = T,delete=F)

# constant survival, constant recapture
Model.9 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.ct, 
                                                    p=p.ct, 
                                                    pent=pent.ct, 
                                                    N=N.time), output = T,delete=F)
# constant survival, time-dependent recapture
Model.10 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.ct,p=p.time, pent=pent.ct, N=N.time),output = T,delete=F)
# time-dependent survival, constant recapture
Model.11 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.time,p=p.ct, pent=pent.ct, N=N.time),output = T,delete=F)
# time-dependent survival, time-dependent recapture
Model.12 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.time,p=p.time, pent=pent.ct, N=N.time),output = T,delete=F)

# constant survival, constant recapture
Model.13 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.ct,p=p.ct, pent=pent.time, N=N.time),output = T,delete=F)
# constant survival, time-dependent recapture
Model.14 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.ct,p=p.time, pent=pent.time, N=N.time),output = T,delete=F)
# time-dependent survival, constant recapture
Model.15 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.time,p=p.ct, pent=pent.time, N=N.time),output = T,delete=F)
# time-dependent survival, time-dependent recapture
Model.16 = mark(hw.proc,hw.ddl,model.parameters=list(Phi=Phi.time,p=p.time, pent=pent.time, N=N.time),output = T,delete=F)


summary(Model.1)$AICc
summary(Model.2)$AICc
summary(Model.3)$AICc
summary(Model.4)$AICc
summary(Model.5)$AICc
summary(Model.6)$AICc
summary(Model.7)$AICc
summary(Model.8)$AICc
# summary(Model.9)$AICc
# summary(Model.10)$AICc
# summary(Model.11)$AICc
# summary(Model.12)$AICc
# summary(Model.13)$AICc
# summary(Model.14)$AICc
# summary(Model.15)$AICc
# summary(Model.16)$AICc


run.caps <- function() {
  
  # Define parameters for p and c 
  # Note: "share=TRUE" indicates that 'p' & 'c'
  #  share the same columns in the design matrix
  Phi.ct = list(formula=~1) # constant
  Phi.time = list(formula=~time) # year effect
  # detection process
  p.ct = list(formula=~1) # constant
  p.time = list(formula=~time) # year effect
  # entryprobablity process
  pent.ct = list(formula=~1)
  pent.time = list(formula=~time)
  # superpopulation process
  # N.ct = list(formula=~1)
  # N.time = list(formula=~time)
  
  # pent.ct = list(formula=~1).   ### doesn't work
  # pent.time = list(formula=~time)### doesn't work
  # 
  # N.ct = list(formuala=~1)### doesn't work
  # N.time = list(formula=~time)### doesn't work
  
  # Create competing models based on strutures for 'p' & 'c'
  caps.model.list = RMark::create.model.list("POPAN")
  
  # NOTE: if you do not want to see the output for each model, add  
  # ', output=FALSE' after 'ddl=caps.ddl' below. 
  caps.results = mark.wrapper(caps.model.list,
                              data = hw.proc, 
                              ddl = hw.ddl)
  
  # Return model table and list of models
  return(caps.results)
}

caps.results <- run.caps()

mark.wrapper(caps.results,
             data = hw.proc, 
             ddl = hw.ddl)

caps.results_df <- as.data.frame(caps.results$model.table)

POPAN_res <- caps.results$Phi.ct.p.ct.pent.ct$results$real
POPAN_res <- rownames_to_column(POPAN_res, "Variable")
write.csv(POPAN_res, "Manuscript_Figures/CJS_Result_tbl.csv", row.names = F)

#write.csv(caps.results_df, "Manuscript_Figures/POPAN_summary_table.csv",  row.names = F)

phitable = get.real(Model.1,"Phi", se= TRUE)
# names(phitable)
phitable[c("estimate","se","lcl","ucl")][1,]

ptable = get.real(Model.1,"p", se= TRUE)
ptable[c("estimate","se","lcl","ucl")][1:11,]

obs = gregexpr("1", ch_df$ch)
n_obs = summary(as.factor(unlist(obs)))
estim_abundance = n_obs[-1]/ptable$estimate[1:10] 
estim_abundance

nb_bootstrap = 100
nb_years = 11
target = data.frame(ch_df,stringsAsFactors=F)
popsize = matrix(NA,nb_bootstrap, nb_years-1)
set.seed(5)
pseudo = target # initialization
#colnames(pseudo)[1] <- "ch"

# define model structure
hw.proc = RMark::process.data(pseudo, model="POPAN")
hw.ddl = RMark::make.design.data(hw.proc)
# define parameter structure
phi.ct = list(formula=~1)
p.time = list(formula=~1)

for (k in 1:nb_bootstrap){
  # resample in the original dataset with replacement
  pseudo$ch = sample(target$ch, replace=T)
  # fit model with Mark
  res = mark(hw.proc,hw.ddl,model.parameters=list(Phi=phi.ct,p=p.time),delete=TRUE,output=FALSE)
  # get recapture prob estimates
  ptable = get.real(res,"p", se= TRUE)
  # calculate the nb of recaptured individiduals / occasion
  allobs = gregexpr("1", pseudo$ch)
  n = summary(as.factor(unlist(allobs)))
  popsize[k,] <- n[-1]/ptable$estimate[1:(nb_years-1)]
}

ci_hw = apply(popsize,2,quantile,probs=c(2.5/100,97.5/100),na.rm=T)
ci_hw

plot(1:(nb_years-1),estim_abundance, col="black", type="n", pch=21, xlab="Years", lty=3, ylab="Estimated abundance", main="dataset 1",lwd=3,ylim=c(0,15))
polygon(c(rev(1:(nb_years-1)), 1:(nb_years-1)), c(rev(ci_hw[2,]), ci_hw[1,]), col = 'grey80', border = NA)
lines(1:(nb_years-1), estim_abundance, col="black",lty=3,type='o',lwd=3,pch=21)

ci_hw_df <- ci_hw %>% 
  t() %>%
  as.data.frame() %>%
  mutate(Year = seq(2015,2024, by=1))

EstAbund <- as.data.frame(estim_abundance) %>%
  #rownames_to_column("Year") %>%
  mutate(Year = seq(2015,2024, by=1)) %>%
  left_join(ci_hw_df)

write.csv(round(EstAbund,3), "Manuscript_Figures/CJS_EstimatedAbundance.csv", row.names = F)

POPAN_plot <- ggplot(EstAbund, aes(x=Year, y=estim_abundance)) + 
  geom_point() +
  geom_line(aes(group=1)) + 
  geom_ribbon(aes(ymin=`2.5%`, 
                              ymax = `97.5%`), alpha=0.3) + 
  scale_x_continuous(breaks=seq(2015,2024,by=1)) + 
  labs( y = "Estimated population size\n(and 95% confidence interval)") + 
  coord_cartesian(ylim=c(-1,15)) + 
  theme(panel.grid.major.y = element_line(linewidth = 0.25))
POPAN_plot

ggsave("Images/STSR presentation/Popan_PopEst_plot.png" , POPAN_plot, 
       width=6, height=4)

#### Calculate net recruitment
Phi  <- get.real(Model.1, "Phi",  se = TRUE)
pent <- get.real(Model.1, "pent", se = TRUE)
Nhat <- Model.1$results$real["N g1 a0 t1", "estimate"]   # row label depends on your model

B <- Nhat * pent$estimate                 # entrants per interval

# Abundance by occasion
N_t    <- numeric(length(B) + 1)
N_t[1] <- Nhat * pent$estimate[1]         # entrants at first occasion (pent0)
for (t in seq_along(Phi$estimate)) {
  N_t[t + 1] <- N_t[t] * Phi$estimate[t] + B[t + 1]
}

# Per-capita recruitment
f <- B[-1] / N_t[-length(N_t)]

