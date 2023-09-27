#version
setwd("C:/Users/shinh/OneDrive/Desktop/nhanes")

library(tidyverse)

# Load Data / Select and Rename Columns 
set_data <- 'nhanes'
folder_path <- './data/'
group_condition <- c('sex','race') #'sex','age_group', 'race'
save_file <- paste0(folder_path,set_data,'_risk_by_sex_race.csv')

ds <- read.csv(paste0(folder_path,set_data,'.csv')) |>
  filter(age < 80) |>
  mutate(age_group = (age %/% 10)  %/%  2 )

str(ds)

# Setting criterion
MALE <- 0
FEMAEL <- 1

GL <- 100
SBP <- 130
DBP <- 85
TG <- 150
HDL_MALE <- 40
HDL_FEMALE <- 50

#KNHANES, NHANES
WC_MALE <- ifelse(set_data=='knhanes',90,102)
WC_FEMALE <- ifelse(set_data=='knhanes',85,88)

ds.new <- mutate(ds,
                 ab_wc = as.integer(ifelse(sex==MALE, 
                                           wc >= WC_MALE, 
                                           wc >= WC_FEMALE)),
                 ab_hdl = as.integer(ifelse(sex==MALE, 
                                            hdl < HDL_MALE, 
                                            hdl < HDL_FEMALE)),
                 ab_tg = as.integer(tg>=TG),
                 ab_gl = as.integer(gl>=GL),
                 ab_bp = as.integer(sbp>=SBP | dbp>=DBP),
                 ab_sum = ab_wc+ab_hdl+ab_tg+ab_gl+ab_bp,
                 mets = as.integer(ab_sum>=3) # mets = 1, non-mets=0
)


# RTAS (Our) 

# step1: Unit Scaling
RATIO = 0.1
unit_scaled_ds <- mutate(ds.new,
                         unit_wc=ifelse(sex==MALE,
                                        (wc-WC_MALE)/(WC_MALE*RATIO),
                                        (wc-WC_FEMALE)/(WC_FEMALE*RATIO)),
                         unit_hdl=ifelse(sex==MALE, 
                                         -1*(hdl-HDL_MALE)/(HDL_MALE*RATIO),
                                         -1*(hdl-HDL_FEMALE)/(HDL_FEMALE*RATIO)),
                         unit_tg=(tg-TG)/(TG*RATIO),
                         unit_gl=(gl-GL)/(GL*RATIO),
                         unit_sbp=(sbp-SBP)/((SBP-DBP)*RATIO),
                         unit_dbp=(dbp-DBP)/((SBP-DBP)*RATIO)
)

str(unit_scaled_ds)

# step2: Apply Sigmoid function

# Elliot Sigmoid Function
elliot <- function(x){
  (0.5*x)/(1+abs(x)) + 0.5
}

sig_ds <- mutate(unit_scaled_ds,
                 sig_wc = elliot(unit_wc),
                 sig_hdl = elliot(unit_hdl),
                 sig_tg = elliot(unit_tg),
                 sig_gl = elliot(unit_gl),
                 sig_sbp = elliot(unit_sbp),
                 sig_dbp = elliot(unit_dbp),
                 sig_bp = ifelse(sig_sbp>sig_dbp,sig_sbp,sig_dbp)
)
head(sig_ds)

# step 3: Generating MetS risk features using three-axis radar charts

# Triangular Area: One Radar Chart
get_triangular_area <- function(x){
  
  p1 <- min(x[1], 0.5)
  p2 <- min(x[2], 0.5)
  p3 <- min(x[3], 0.5)
  
  s_in <- 0.5*sin(pi*1/3)*(p1*p2+p2*p3+p3*p1)
  s_out <- 0.5*sin(pi*1/3)*(x[1]*x[2]+x[2]*x[3]+x[3]*x[1])
  s_max <- 0.5*sin(pi*1/3)*(0.75) # x1=x2=x3=0.5
  
  s1 <- s_in/s_max
  s2 <- s_out/s_max
  s1_max <- s_max/s_max
  
  s_alpha <- round(s1_max - s1,8)
  s_beta <- s2-s1
  s_left <- (s1+s_beta)/(1+s_beta)
  s_right <- s_beta/3
  
  i <- ifelse(s_alpha == 0, 1, 0)
  s3 <- 0.5*(s_left + i*s_right)
  
  return(s3)
}

# step 4: MetS Risk Score : Ten Radar Chart
get_risk_score <- function(sig_ds){
  
  chart_axis <- combn(c('sig_wc','sig_gl','sig_tg','sig_hdl','sig_bp'),3) 
  
  score <- c()
  for(i in 1:ncol(chart_axis)){
    s3 <- sig_ds |> select(chart_axis[,i]) |> apply(1,get_triangular_area)
    score <- cbind(score,s3)
  }
  risk_score <- sqrt(apply(score,1,mean))
  risk_type <- apply(score,1,function(x){which(x==max(x))[1]}) # First choice if multiple exist
  
  # Type info: 
  # 1 = wc gl tg, 2 = wc gl hdl, 3 = wc gl bp, 4 = wc tg hdl, 5 = wc tg bp
  # 6 = wc hdl bp, 7 = gl tg hdl, 8 = gl tg bp, 9 = gl hdl bp, 10 = tg hdl bp
  
  return(cbind(RTAS=risk_score,RTAS_type=risk_type))
}

risk_type <- get_risk_score(sig_ds)

rtas <- cbind(sig_ds,risk_type)
str(rtas)



# siMS score

cols <- c('id','sex','race','age_group','height','wc','gl','tg','sbp','hdl','ab_sum','mets')

siMS <- ds.new |>
  select(all_of(cols)) |>
  mutate(siMS = 2*wc/height + 
           gl/GL + 
           sbp/SBP - 
           hdl/ifelse(sex==0,HDL_MALE,HDL_FEMALE))

head(siMS)



# cMets
# step 1: MAP, (SBP-DBP)/3 + DBP
# step 2: Standardization
# step 3: Regression, Each factor ~ Age + Sex or Age + Sex + Race
# step 4: Residual, Standardization

#group_condition <- c('sex') #,'age_group','race'

scale_this <- function(x) as.vector(scale(x))

cMetS <- ds.new |>
  # step 1: MAP, (SBP-DBP)/3 + DBP
  mutate(map = (sbp-dbp)/3+dbp) |> 
  group_by(across(all_of(group_condition))) |>
  # step 2: Standardization
  mutate(across(c('map','wc','gl','tg','hdl'), scale_this))

# step 3: Regression, Each factor ~ Age + Sex or Age + Sex + Race
map_rsd <- cMetS |> group_split() |> map(~lm(map~age, data=.)) |> map('residuals') 
wc_rsd <- cMetS |> group_split() |> map(~lm(wc~age, data=.)) |> map('residuals') 
gl_rsd <- cMetS |> group_split() |> map(~lm(gl~age, data=.)) |> map('residuals') 
tg_rsd <- cMetS |> group_split() |> map(~lm(tg~age, data=.)) |> map('residuals') 
hdl_rsd <- cMetS |> group_split() |> map(~lm(hdl~age, data=.)) |> map('residuals') 
  
# step 4: Residual, Standardization
ds_all <- c()
ds_rsd <- c()
group_info <- attributes(cMetS)$group
str(cMetS)
n <- ifelse(is.null(nrow(group_info)),1,nrow(group_info))
for(i in 1:n){
  if(n != 1){
    ds_grp <- cMetS[group_info$.rows[[i]],]  
  }else{
    ds_grp <- cMetS
  }
  ds_rsd <- cbind(map_rsd=map_rsd[[i]], 
                  wc_rsd=wc_rsd[[i]], 
                  gl_rsd=gl_rsd[[i]], 
                  tg_rsd=tg_rsd[[i]], 
                  hdl_rsd=hdl_rsd[[i]])
  ds_rsd <- cbind(ds_grp, ds_rsd)
  ds_all <- rbind(ds_all, ds_rsd) |> as_tibble()
  ds_rsd <- c()
}

str(ds_all)

cMetS <- ds_all |> 
  group_by(across(all_of(group_condition))) |>
  mutate(across(c('map_rsd','wc_rsd','gl_rsd','tg_rsd','hdl_rsd'), scale_this)) |>
  # step 5: Sum Residual
  mutate(cMetS = map_rsd+wc_rsd+gl_rsd+tg_rsd-hdl_rsd) |>
  ungroup()

str(cMetS)




# MetSSS
# https://cran.r-project.org/web/packages/pscore/vignettes/metsss.html
#library(pscore)

apply_censor <- function(x, y) ifelse(x <= y, 0, x-y)
apply_rev_censor <- function(x, y) ifelse(x > y, 0, y-x)

scale_this <- function(x) as.vector(scale(x))

MetSSS <- ds.new |>
  # step 1: censor & center
  mutate(wc = apply_censor(wc, ifelse(sex==0, WC_MALE, WC_FEMALE)),
         sbp = apply_censor(sbp, SBP),
         dbp = apply_censor(dbp, DBP),
         tg = apply_censor(tg, TG),
         gl = apply_censor(gl, GL),
         hdl = apply_rev_censor(hdl, ifelse(sex==0, HDL_MALE, HDL_FEMALE)),
  ) |> 
  group_by(across(all_of(group_condition))) |>
  # step 2: Standardization
  mutate(across(c('wc','sbp','dbp','gl','tg','hdl'), scale_this)) 

# step 3: Orthogonalize
pca_matirx <- MetSSS |>
  group_split() |>  
  map(~prcomp(.|> select(wc,sbp,dbp,gl,tg,hdl))) |>
  map('rotation')

sdev <- MetSSS |>
  group_split() |>  
  map(~prcomp(.|> select(wc,sbp,dbp,gl,tg,hdl))) |>
  map('sdev')

ds_all <- c()
ds_rsd <- c()
group_info <- attributes(MetSSS)$group

n <- ifelse(is.null(nrow(group_info)),1,nrow(group_info))

for(i in 1:n){
  if(n != 1){
    ds_grp <- MetSSS[group_info$.rows[[i]],] |> ungroup()
  }else{
    ds_grp <- MetSSS |> ungroup()
  }
  pca_grp <- pca_matirx[[i]]
  
  s2 <- diag(1/sdev[[i]])
  l <- pca_grp
  
  x <- ds_grp |> select(wc,sbp,dbp,gl,tg,hdl) |> as.matrix()  
  
  
  # step 4: standardize
  res <- x %*% (l %*% s2)
  
  # step 5: square, sum, square root
  risk_score = res^2 |> rowMeans() |> sqrt()
  
  ds_rsd <- cbind(ds_grp, MetSSS=risk_score)
  ds_all <- rbind(ds_all, ds_rsd)
  ds_rsd <- c()
}

MetSSS <- ds_all
str(MetSSS)



# ASD
min_max_scale <- function(x){(x-min(x))/(max(x)-min(x))}
threshold_scale <- function(x,threshold){(threshold-min(x))/(max(x)-min(x))}

ds.new <- ds |>
  group_by(across(all_of(group_condition))) |>
  # step 1: Normalize(min-max scaling), 
  mutate(wc_s = min_max_scale(wc),
         sbp_s = min_max_scale(sbp),
         dbp_s = min_max_scale(dbp),
         tg_s = min_max_scale(tg),
         gl_s = min_max_scale(gl),
         hdl_s = 1-min_max_scale(hdl),
         hdl_min = min(hdl),
         hdl_max = max(hdl)
  ) |>
  # scaled thresholds
  mutate(wc_th = threshold_scale(wc,ifelse(sex==MALE,WC_MALE,WC_FEMALE)),
         hdl_th = 1-threshold_scale(hdl,ifelse(sex==MALE,HDL_MALE,HDL_FEMALE)),
         sbp_th = threshold_scale(sbp,SBP),
         dbp_th = threshold_scale(dbp,DBP),
         tg_th = threshold_scale(tg,TG),
         gl_th = threshold_scale(gl,GL),
  ) |>
  mutate(bp_s = ifelse(sbp_s > dbp_s, sbp_s, dbp_s),
         bp_th = ifelse(sbp_s==bp_s, sbp_th, dbp_th)
  ) |>
  mutate(ab_wc = as.integer(wc_s>=wc_th),
         ab_hdl = as.integer(hdl_s>hdl_th),
         ab_tg = as.integer(tg_s>=tg_th),
         ab_gl = as.integer(gl_s>=gl_th),
         ab_sbp = as.integer(sbp_s>=sbp_th),
         ab_dbp = as.integer(dbp_s>=dbp_th),
         ab_bp = as.integer(bp_s>=bp_th),
         ab_sum = ab_wc+ab_hdl+ab_tg+ab_gl+ifelse(ab_sbp+ab_dbp != 0, 1, 0),
         mets = as.integer(ab_sum>=3) # mets = 1, non-mets=0
  )

# step 2: weight
# Number of thresholds exceeded in one risk factor / Sum of the number of thresholds exceeded for all risk factors

weight <- ds.new |> 
#  group_by(across(all_of(group_condition))) |>
  do(count= colSums(.|> select(ab_wc,ab_sbp,ab_dbp,ab_tg,ab_gl,ab_hdl)))

weights <- c()
group_info <- attributes(ds.new)$groups

n <- ifelse(is.null(nrow(group_info)),1,length(weight$count))

for(i in 1:n){
  w <- weight$count[[i]]
  ab_bp <- ifelse(w['ab_sbp']>w['ab_dbp'],w['ab_sbp'],w['ab_dbp']) # The larger one is selected
  w <- c(w[c('ab_wc','ab_tg','ab_gl','ab_hdl')],ab_bp)
  names(w) <- c('wc_w','tg_w','gl_w','hdl_w','bp_w')
  w <- c(unlist(group_info[i,-length(group_info)]),w/sum(w))
  
  weights <- rbind(weights,w)
}
weights <- as_tibble(weights)
str(weights)
# step 3: quantify function

get_intersection_area <- function(x1,x2,t1,t2){
  s<-0
  if(x1>=t1 & x2>=t2){
    s <- 1
  }else if(x1<t1 & x2<t2){
    s <- x1*x2/t1*t2
  }else if(x1>t1 & x2<t2){
    q <- x1*t1*(t2-x2)^2/(x2*(x1-t1)+x1*(t2-x2))
    s <- 1-q/t1*t2
  }else if(x1<t1 & x2>t2){
    q <- x2*t2*(t1-x1)^2/(x1*(x2-t2)+x2*(t1-x1))
    s <- 1-q/t1*t2
  }
  return(s)
}

# step 4: weighted sum, five polygon
if(n != 1){
  asd <- ds.new |> ungroup() |>
    left_join(weights, by=group_condition)  
}else{
  w <- matrix(rep(unlist(weights), times = nrow(ds)), 
         byrow = T, 
         nrow = nrow(ds),
         dimnames = list(NULL,
                         c('wc_w','tg_w','gl_w','hdl_w','bp_w'))
  )
  asd <- cbind(ds.new,w)
}

#the order of axis: wc, gl, bp, tg, hdl
risk_score <- c()
for(i in 1:nrow(asd)){
  ins <- asd[i,]
  wc_risk <- get_intersection_area(ins$wc_s,ins$hdl_s,ins$wc_th,ins$hdl_th)
  gl_risk <- get_intersection_area(ins$gl_s,ins$wc_s,ins$gl_th,ins$wc_th)
  bp_risk <- get_intersection_area(ins$bp_s,ins$gl_s,ins$bp_th,ins$gl_th)
  tg_risk <- get_intersection_area(ins$tg_s,ins$bp_s,ins$tg_th,ins$bp_th)
  hdl_risk <-get_intersection_area(ins$hdl_s,ins$tg_s,ins$hdl_th,ins$tg_th)
  
  risk_score[i]=wc_risk*ins$wc_w + gl_risk*ins$gl_w + bp_risk*ins$bp_w + 
    tg_risk*ins$tg_w + hdl_risk*ins$hdl_w
}

asd <- cbind(asd,ASD=risk_score)
str(asd)



### ----------------------------------------------
### all risk scores
### ----------------------------------------------
head(ds)

ds.new <- mutate(ds,
                 ab_wc = as.integer(ifelse(sex==MALE, 
                                           wc >= WC_MALE, 
                                           wc >= WC_FEMALE)),
                 ab_hdl = as.integer(ifelse(sex==MALE, 
                                            hdl < HDL_MALE, 
                                            hdl < HDL_FEMALE)),
                 ab_tg = as.integer(tg>=TG),
                 ab_gl = as.integer(gl>=GL),
                 ab_bp = as.integer(sbp>=SBP | dbp>=DBP),
                 ab_sum = ab_wc+ab_hdl+ab_tg+ab_gl+ab_bp,
                 mets = as.integer(ab_sum>=3) # mets = 1, non-mets=0
)

risk_score <- ds.new |> 
  left_join(siMS |> select(id, siMS), by='id') |>
  left_join(cMetS |> select(id, cMetS), by='id') |>
  left_join(MetSSS |> select(id, MetSSS), by='id') |>
  left_join(asd |> select(id, ASD), by='id') |>
  left_join(rtas |> select(id, RTAS, RTAS_type), by='id') |>
  select(id,sex,race,age_group,age,
         wc,height,sbp,dbp,gl,tg,hdl,
         ab_wc,ab_hdl,ab_tg,ab_gl,ab_bp,
         siMS,cMetS,MetSSS,ASD,RTAS,RTAS_type,
         ab_sum,mets)

head(risk_score)
write.csv(risk_score,save_file,row.names = F)
