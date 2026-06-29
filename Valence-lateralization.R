###Emotional Valence Lateralisation 27/June 2026

##
rm(list=ls())
library(boot)
library(lme4);library(lmerTest)
library(car)
library(DHARMa)
library(emmeans)
sem <- function(x){sd(x,na.rm=T)/sqrt(length(x[!is.na(x)]))} 

data <- read.csv("emval-data.csv",stringsAsFactors=FALSE)
attach(data)
names(data)
#[1] "fish.id"   "sex"       "stimulus"  "st.object" "order"     "weight"    "sl"        "right"     "left"      "LI"     

abs.LI <- abs(LI)
dat.new <- cbind(data,abs.LI)
rm(abs.LI)
detach(data);attach(dat.new)
length(unique(fish.id))
#[1] 65

### excluding LI = ±1 data points
dat.new.c <- subset(dat.new, abs(LI)!=1)
detach(dat.new);attach(dat.new.c)
length(unique(fish.id))
#[1] 61

###################################


# Population level Lateralisation

t.test(LI[stimulus=="pos" & sex=="F"])
#t = 0.84589, df = 24, p-value = 0.406
t.test(LI[stimulus=="neg" & sex=="F"])
#t = 0.2997, df = 24, p-value = 0.767
t.test(LI[stimulus=="pos" & sex=="M"])
#t = -0.51402, df = 30, p-value = 0.611
t.test(LI[stimulus=="neg" & sex=="M"])
#t = 0.97276, df = 24, p-value = 0.3404

mean.b <- function(data.in,i) {mean(data.in$LI[i],na.rm=T)}
boot.ci(boot(dat.new.c[stimulus=="pos" & sex=="F",],mean.b,R=1000),type="bca") #intervals calculated using the adjusted bootstrap percentile (BCa) method.
boot.ci(boot(dat.new.c[stimulus=="neg" & sex=="F",],mean.b,R=1000),type="bca")
boot.ci(boot(dat.new.c[stimulus=="pos" & sex=="M",],mean.b,R=1000),type="bca")
boot.ci(boot(dat.new.c[stimulus=="neg" & sex=="M",],mean.b,R=1000),type="bca")

# Individual level lateralisation

t.test(abs.LI[stimulus=="pos" & sex=="F"])
#t = 5.4701, df = 24, p-value = 1.271e-05
t.test(abs.LI[stimulus=="neg" & sex=="F"])
#t = 6.5133, df = 24, p-value = 9.78e-07
t.test(abs.LI[stimulus=="pos" & sex=="M"])
#t = 8.9766, df = 30, p-value = 5.312e-10
t.test(abs.LI[stimulus=="neg" & sex=="M"])
#t = 5.6461, df = 24, p-value = 8.184e-06

mean.b <- function(data.in,i) {mean(data.in$abs.LI[i],na.rm=T)}
boot.ci(boot(dat.new.c[stimulus=="pos" & sex=="F",],mean.b,R=1000),type="bca") #intervals calculated using the adjusted bootstrap percentile (BCa) method.
boot.ci(boot(dat.new.c[stimulus=="neg" & sex=="F",],mean.b,R=1000),type="bca") #intervals
boot.ci(boot(dat.new.c[stimulus=="pos" & sex=="M",],mean.b,R=1000),type="bca") #intervals
boot.ci(boot(dat.new.c[stimulus=="neg" & sex=="M",],mean.b,R=1000),type="bca") #intervals

#lmer

model_LI <- lmer(LI ~ stimulus*sex + (1|fish.id))
resids <- simulateResiduals(model_LI)

# Plot diagnostics
# This gives you a QQ-plot (Normality) and Residual vs. Predicted (Homoscedasticity)
plot(resids)

# Specifically test for outliers and heteroscedasticity (needs both p > .05 for non-violation of assumptions)
testOutliers(resids)
testDispersion(resids)

summary(model_LI) # Assumptions are met
#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']

#Fixed effects:
#                 Estimate.  Std. Error  df     t value Pr(>|t|)
#(Intercept)       0.02759    0.08781 93.34645   0.314    0.754
#stimuluspos       0.01192    0.09939 46.25108   0.120    0.905
#sexM              0.01255    0.12351 95.96279   0.102    0.919
#stimuluspos:sexM -0.10346    0.13788 47.48730  -0.750    0.457

#######

model <- lmer((abs.LI) ~ stimulus*sex + (1|fish.id))

resids <- simulateResiduals(model)

# Plot diagnostics
# This gives you a QQ-plot (Normality) and Residual vs. Predicted (Homoscedasticity)
plot(resids)

# Specifically test for outliers and heteroscedasticity (needs both p > .05 for non-violation of assumptions)
testOutliers(resids)
testDispersion(resids)

summary(model) # Assumptions are met
#Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']

#Fixed effects:
#                 Estimate Std. Error     df    t value Pr(>|t|)    
#(Intercept)       0.35426    0.05194 97.39417   6.821 7.68e-10 ***
#stimuluspos      -0.09904    0.06270 49.68991  -1.580  0.12056    
#sexM             -0.06083    0.07322 98.92867  -0.831  0.40816    
#stimuluspos:sexM  0.25676    0.08687 51.06905   2.956  0.00471 ** 

###########

#Simple effects analysis

## comparing levels of stimulus within each levels of sex
stim.emm<- emmeans(model, ~ stimulus | sex) 
                        
pairs(stim.emm)

#sex = F:
#  contrast  estimate     SE   df t.ratio p.value
#neg - pos    0.099 0.0629 50.5   1.575  0.1216

#sex = M:
#  contrast  estimate     SE   df t.ratio p.value
# neg - pos   -0.158 0.0604 53.4  -2.612  0.0117

#Degrees-of-freedom method: kenward-roger

#plot(stim.emm, comparisons = TRUE)


eff_size(stim.emm, sigma(model), df.residual(model))

#sex = F:
#contrast  effect.size    SE   df lower.CL upper.CL
#neg - pos       0.455 0.291 97.5   -0.122    1.032

#sex = M:
#  contrast  effect.size    SE   df lower.CL upper.CL
#neg - pos      -0.725 0.282 96.3   -1.285   -0.165



## comparing levels of sex within each levels of stimulus
sex.emm <- emmeans(model, pairwise ~ sex | stimulus)
                        
pairs(sex.emm)
#stimulus = neg:
#contrast estimate     SE df t.ratio p.value
#F - M      0.0608 0.0735 99   0.828  0.4097

#stimulus = pos:
#  contrast estimate     SE df t.ratio p.value
#F - M     -0.1959 0.0700 97  -2.799  0.0062

#Degrees-of-freedom method: kenward-roger 

#plot(sex.emm, comparisons = TRUE)
## 

eff_size(sex.emm, sigma(model), df.residual(model))
#stimulus = neg:
#contrast effect.size    SE   df lower.CL upper.CL
#F - M           0.28 0.338 97.5   -0.392    0.951

#stimulus = pos:
#  contrast effect.size    SE   df lower.CL upper.CL
#F - M          -0.90 0.328 96.3   -1.551   -0.250

## 
# Do we need to correct for multiple tests here?

#combined <- rbind(pairs(stim.emm), pairs(sex.emm))
#summary(combined, adjust = "holm")


emm <- emmeans(model, ~ stimulus * sex)
contrast <- pairs(emm, adjust = "tukey")
contrast
#contrast      estimate     SE    df t.ratio p.value
#neg F - pos F   0.0990 0.0629 50.5   1.575  0.4020
#neg F - neg M   0.0608 0.0735 99.0   0.828  0.8411
#neg F - pos M  -0.0969 0.0700 97.0  -1.384  0.5121
#pos F - neg M  -0.0382 0.0735 99.0  -0.520  0.9541
#pos F - pos M  -0.1959 0.0700 97.0  -2.799  0.0309
#neg M - pos M  -0.1577 0.0604 53.4  -2.612  0.0550

#Degrees-of-freedom method: kenward-roger 
#P value adjustment: tukey method for comparing a family of 4 estimates  


## Graphs

means <- c(mean(LI[stimulus=="pos" & sex=="F"],na.rm=T),
           mean(LI[stimulus=="pos" & sex=="M"],na.rm=T),
           mean(LI[stimulus=="neg" & sex=="F"],na.rm=T),
           mean(LI[stimulus=="neg" & sex=="M"],na.rm=T))

sems <- c(sem(LI[stimulus=="pos" & sex=="F"]),
          sem(LI[stimulus=="pos" & sex=="M"]),
          sem(LI[stimulus=="neg" & sex=="F"]),
          sem(LI[stimulus=="neg" & sex=="M"]))
bars <- barplot(means,
                ylim=c(-1,1),
                col=c("white","grey30","white","grey30"),
                names.arg=c("Female","Males","Female","Males"),
                ylab="Mean Laterality Index", 
                #main = "Laterality (excluding LI= ±1) for Emotional Valence"
)

arrows(bars,means-sems,bars,means+sems,
       length = 0.10, # width of the arrowhead
       angle = 90, # angle of the arrowhead
       code = 3 ## arrowhead in both ends
       )


# add the bootstrap 95% CI
lines(c(bars[1]+0.25,bars[1]+0.25),c(-0.0649, 0.1984),lwd=5,col="grey")
lines(c(bars[2]+0.25,bars[2]+0.25),c(-0.1337, 0.2178),lwd=5,col="grey")
lines(c(bars[3]+0.25,bars[3]+0.25),c(-0.2396, 0.1281),lwd=5,col="grey")
lines(c(bars[4]+0.25,bars[4]+0.25),c(-0.2396, 0.1281),lwd=5,col="grey")

# Now, to add the individual data points with jitter:

grp <- list(
  LI[stimulus=="pos" & sex=="F"],
  LI[stimulus=="pos" & sex=="M"],
  LI[stimulus=="neg" & sex=="F"],
  LI[stimulus=="neg" & sex=="M"]
)

for (i in 1:4) {
  points(
    jitter(rep(bars[i], length(grp[[i]])), amount = 0.15), 
    grp[[i]],                                              
    pch = 16,                                                  
    col = adjustcolor(c("snow4", "black", "snow4", "black")[i] , alpha.f = 0.8
                      ),
    cex = 0.8                                                   
  )
}


lines(c(0,7),c(0,0))
mtext("Positive Valence",side=1,line=3,at=(bars[1]+bars[2])/2, cex=1)
mtext("Negative Valence",side=1,line=3,at=(bars[3]+bars[4])/2, cex=1)


text(bars[1],-1,paste("N=",sum(!is.na(LI[stimulus=="pos" & sex=="F"]))),pos=3)
text(bars[2],-1,paste("N=",sum(!is.na(LI[stimulus=="pos" & sex=="M"]))),pos=3)
text(bars[3],-1,paste("N=",sum(!is.na(LI[stimulus=="neg" & sex=="F"]))),pos=3)
text(bars[4],-1,paste("N=",sum(!is.na(LI[stimulus=="neg" & sex=="M"]))),pos=3)
dev.print(device=pdf,file="Fig1-LI.pdf")


#ALI

means <- c(mean(abs.LI[stimulus=="pos" & sex=="F"],na.rm=T),
           mean(abs.LI[stimulus=="pos" & sex=="M"],na.rm=T),
           mean(abs.LI[stimulus=="neg" & sex=="F"],na.rm=T),
           mean(abs.LI[stimulus=="neg" & sex=="M"],na.rm=T))

sems <- c(sem(abs.LI[stimulus=="pos" & sex=="F"]),
          sem(abs.LI[stimulus=="pos" & sex=="M"]),
          sem(abs.LI[stimulus=="neg" & sex=="F"]),
          sem(abs.LI[stimulus=="neg" & sex=="M"]))


bars <- barplot(means,
                
                ylim=c(0,1),
                col=c("white","grey30","white","grey30"),
                names.arg=c("Female","Males","Female","Males"),
                ylab="Absolute Laterality Index",
                
                cex.names=1,
                cex.axis = 1.2,         # Increase size of axis tick labels
                font.axis = 1.5)       # Increases size of x-axis category names

arrows(bars,means-sems,bars,means+sems,
       length = 0.10, # width of the arrowhead
       angle = 90, # angle of the arrowhead
       code = 3 ## arrowhead in both ends
)

# add the bootstrap 95% CI
lines(c(bars[1]+0.25,bars[1]+0.25),c(0.1799, 0.3659),lwd=5,col="grey")
lines(c(bars[2]+0.25,bars[2]+0.25),c(0.2525, 0.4684),lwd=5,col="grey")
lines(c(bars[3]+0.25,bars[3]+0.25),c(0.3585, 0.5585),lwd=5,col="grey")
lines(c(bars[4]+0.25,bars[4]+0.25),c(0.1901, 0.3851),lwd=5,col="grey")


# Now, to add the individual data points without jitter:
grp <- list(
  abs.LI[stimulus=="pos" & sex=="F"],
  abs.LI[stimulus=="pos" & sex=="M"],
  abs.LI[stimulus=="neg" & sex=="F"],
  abs.LI[stimulus=="neg" & sex=="M"]
)

for (i in 1:4) {
  points(
    jitter(rep(bars[i], length(grp[[i]])), amount = 0.15), 
    grp[[i]],                                              
    pch = 16,                                                  
    col = adjustcolor(c("snow4", "black", "snow4", "black")[i] , #alpha.f = 0.8
    ),
    cex = 0.8                                                   
  )
}

lines(c(0,7),c(0,0))
mtext("Positive Valence",side=1,line=3,at=(bars[1]+bars[2])/2, cex=1)
mtext("Negative Valence",side=1,line=3,at=(bars[3]+bars[4])/2, cex=1)

lines(c(bars[2],bars[4]),c(0.65,0.65), lwd = 2)
lines(c(bars[2],bars[1]),c(0.55,0.55), lwd = 2)

text(((bars[4]+bars[2])/2),0.67,"*",cex=2.5)
text(((bars[1]+bars[2])/2),0.57,"**",cex=2.5)

text(bars[1],0,paste("N=",sum(!is.na(abs.LI[stimulus=="pos" & sex=="F"]))),pos=3)
text(bars[2],0,paste("N=",sum(!is.na(abs.LI[stimulus=="pos" & sex=="M"]))),pos=3)
text(bars[3],0,paste("N=",sum(!is.na(abs.LI[stimulus=="neg" & sex=="F"]))),pos=3)
text(bars[4],0,paste("N=",sum(!is.na(abs.LI[stimulus=="neg" & sex=="M"]))),pos=3)
dev.print(device=pdf,file="Fig2-ALI.pdf")




