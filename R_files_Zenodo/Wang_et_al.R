################################################################################
########           Wang_et_al._2020_R_Codes_for_Calculations            ########
################################################################################

# quartz(width=5, height=7); par(mfrow=c(3, 2))

rm(list=ls())
f_path = ""
setwd(f_path)

com1 = read.csv("Wang_et_al._2020_Raw_Community_Dataset.csv", header=T, fileEncoding="GBK")
cli1 = read.csv("Wang_et_al._2020_Raw_Climate_Dataset.csv", header=T, fileEncoding="GBK")

################################################################################
# The following codes were used to calculate targeted variables, such as 
# multiple-year mean growing-season precipitation, community temporal CV, 
# weighted average species temporal CV and species synchrony and so on.
################################################################################

################################################################################
# Following codes used to calculate:
# bm_site: multiple-year mean community biomass of each site
# pre_site: mean growing-season precipitation of each site
# pre_cv: CV of growing-season precipitationtemperature of each site

com1$cBM = ifelse(apply(com1[, 7:251], 1, sum, na.rm=T)==0, NA, apply(com1[, 7:251], 1, sum, na.rm=T))
biomass_label = ifelse(is.na(tapply(com1$cBM, list(com1$Site, com1$Year), mean, na.rm=T)), NA, 1)

pre = tapply(cli1$Prep, list(cli1$Site, cli1$Year), sum, na.rm=T)
for (i in 1:5) {
    pre[, i] = pre[, i] * biomass_label[, i]
}
pre_site = apply(pre, 1, function(x)(mean(x, na.rm=T)))
pre_cv = apply(pre, 1, function(x)(sd(x, na.rm=T)))/pre_site
bm_site = tapply(com1$cBM, list(com1$Site), mean, na.rm=T)

################################################################################
# Following codes used to calculate:
# sync_site: species synchrony of each site

spe_bm_year = NULL
for (i in 1:length(unique(com1$Site))) {
    temp1 = subset(com1, Site==unique(com1$Site)[i])
    
    for (j in 1:length(unique(temp1$Year))) {
        temp2 = subset(temp1, Year==unique(temp1$Year)[j])
        
        spe_bm_year1 = temp2[1, ]
        spe_bm_year1[1, 7:251] = apply(temp2[, 7:251], 2, function(x)(sum(x, na.rm=T)))/3
        spe_bm_year1[1, ] = ifelse(spe_bm_year1[1,]==0, NA, spe_bm_year1[1,])
        
        spe_bm_year = rbind(spe_bm_year, spe_bm_year1)
    }
}

spe_bmsd_site = as.data.frame(matrix(NA, ncol=251, nrow=23, dimnames =list(NULL, c(names(spe_bm_year)[1:251]))))
for (i in 1:23) {
    spe_bmsd_site[i, 1:6] = spe_bm_year[i, 1:6]
    for (j in 7:251) {
        spe_bmsd_site[i, j] = sd(spe_bm_year[spe_bm_year$Site==i, j], na.rm=T)
    }
}

bmvar_site = apply(tapply(com1$cBM, list(com1$Site, com1$Year), mean, na.rm=T), 1, function(x)var(x, na.rm=T))
sumvar_site = (apply(spe_bmsd_site[, 7:251], 1, function(x)sum(x, na.rm=T))^2)
sync_site = (bmvar_site / ifelse(sumvar_site==0, NA, sumvar_site))

################################################################################
# Following codes used to calculate:
# cbm_cv: community temporal CV of each site

cbm_cv = sqrt(bmvar_site)/bm_site

################################################################################
# Following codes used to calculate:
# spe_rsb: relative species biomass of individual species of each site
# spe_bm: mean species biomass individual species of each site

spe_rsb = NULL
spe_bm = NULL
for (i in 1:23) {
    # i = 1
    temp1 = subset(com1, Site==unique(Site)[i])
    
    spe_rsb1 = temp1[1, ]
    spe_bm1 = temp1[1, ]
    
    total_bm = sum(temp1[, 252], na.rm=T)
    
    if (total_bm!=0) {
        spe_rsb1[, 7:252] = 100 * apply(temp1[, 7:252], 2, function(x)sum(x, na.rm=T))/total_bm
        spe_bm1[, 7:252] = apply(temp1[, 7:252], 2, function(x)mean(x, na.rm=T))
    } else {
        spe_rsb1[1, ] = NA
        spe_bm1[1, ] = NA
    }
    
    spe_rsb1[1, ] = ifelse(spe_rsb1[1, ]==0, NA, spe_rsb1[1, ])
    spe_bm1[1, ] = ifelse(spe_bm1[1, ]==0, NA, spe_bm1[1, ])
    spe_rsb = rbind(spe_rsb, spe_rsb1)
    spe_bm = rbind(spe_bm, spe_bm1)
}

################################################################################
# Following codes used to calculate:
# w_spe_cv: temporal CVs of indivicual species of each site that weighted by 
#           their contribution to community (site) total community biomass
# spe_var: temporal variance of individual species of each site
################################################################################
# Then, used to calculate the weighted average species temporal CV and its components
# of different abundance groups
################################################################################

w_spe_cv = NULL
spe_var = NULL
for (i in 1:23) {
    # i = 1
    temp1 = subset(spe_bm_year, Site==unique(Site)[i])
    
    w_spe_cv1 = temp1[1, ]
    spe_var1 = temp1[1, ]
    mean_cbm = mean(temp1[, 252], na.rm=T)
    
    if (mean_cbm!=0 & !is.na(mean_cbm)) {
        w_spe_cv1[, 7:251] = apply(temp1[, 7:251], 2, function(x)sd(x, na.rm=T))/mean_cbm
        spe_var1[, 7:251] = apply(temp1[, 7:251], 2, function(x)var(x, na.rm=T))
    } else {
        w_spe_cv1[1, ] = NA
        spe_var1[1, ] = NA
    }
    
    w_spe_cv = rbind(w_spe_cv, w_spe_cv1)
    spe_var = rbind(spe_var, spe_var1)
}

################################################################################
# Following codes used to calculate:
# w_specv_site: weighted average species temporal CV of each site

w_specv_site = apply(w_spe_cv[, 7:251], 1, function(x)sum(x, na.rm=T))

################################################################################
# Following codes used to calculate:
# w_spedom_cv: weighted average dominant species temporal CV of each site
# dom_spesd_site: standard deviation of dominant species group of each site
# dom_spebm_site: biomass of dominant species group of each site
# w_specom_cv: weighted average common species temporal CV of each site
# w_sperare_cv: weighted average rare species temporal CV of each site

w_spedomcv_site1 = w_spe_cv
dom_spesd_site1 = spe_var
dom_spebm_site1 = spe_bm
spedom_bool = spe_rsb
for (i in 1:23) {
    spedom_bool[i, 7:251] = ifelse(spe_rsb[i, 7:251]>=5, 1, NA)
}

for (i in 1:23) {
    temp_spedom_bool1 = spedom_bool[i, 7:251]
    temp_w_spe_cv1 = w_spe_cv[i, 7:251]
    temp_w_spe_var1 = spe_var[i, 7:251]
    temp_dombm1 = spe_bm[i, 7:251]
    
    temp_spedom_bool = as.matrix(temp_spedom_bool1)
    temp_w_spe_cv = as.matrix(temp_w_spe_cv1)
    temp_w_spe_var = as.matrix(temp_w_spe_var1)
    temp_dombm = as.matrix(temp_dombm1)
    
    w_spedomcv_site1[i, 7:251] = temp_w_spe_cv * temp_spedom_bool
    dom_spesd_site1[i, 7:251] = sqrt(temp_w_spe_var) * temp_spedom_bool
    dom_spebm_site1[i, 7:251] = temp_dombm * temp_spedom_bool
}

w_spedomcv_site = apply(w_spedomcv_site1[, 7:251], 1, function(x)sum(x, na.rm=T))
dom_spesd_site = apply(dom_spesd_site1[, 7:251], 1, function(x)sum(x, na.rm=T))
dom_spebm_site = apply(dom_spebm_site1[, 7:251], 1, function(x)sum(x, na.rm=T))

########

w_specomcv_site1 = w_spe_cv
specom_bool = spe_rsb
for (i in 1:23) {
    specom_bool[i, 7:251] = ifelse(spe_rsb[i, 7:251]<5 & spe_rsb[i, 7:251]>1, 1, NA)
}

for (i in 1:23) {
    # i = 1
    temp_specom_bool1 = specom_bool[i, 7:251]
    temp_w_spe_cv1 = w_spe_cv[i, 7:251]
    
    temp_specom_bool = as.matrix(temp_specom_bool1)
    temp_w_spe_cv = as.matrix(temp_w_spe_cv1)

    w_specomcv_site1[i, 7:251] = temp_w_spe_cv * temp_specom_bool
}

w_specomcv_site = apply(w_specomcv_site1[, 7:251], 1, function(x)sum(x, na.rm=T))

########

w_sperarecv_site1 = w_spe_cv
sperare_bool = spe_rsb
for (i in 1:23) {
    sperare_bool[i, 7:251] = ifelse(spe_rsb[i, 7:251]<1, 1, NA)
}

for (i in 1:23) {
    temp_sperare_bool1 = sperare_bool[i, 7:251]
    temp_w_spe_cv1 = w_spe_cv[i, 7:251]
    temp_sperare_bool = as.matrix(temp_sperare_bool1)
    temp_w_spe_cv = as.matrix(temp_w_spe_cv1)
    
    w_sperarecv_site1[i, 7:251] = temp_w_spe_cv * temp_sperare_bool
}

w_sperarecv_site = apply(w_sperarecv_site1[, 7:251], 1, function(x)sum(x, na.rm=T))

################################################################################
# Following codes used to calculate:
# richCMUN: community species richness of each site
# richDOMI: dominant species richness of each site
# richCOMM: common species richness of each site
# richRARE: rare species richness of each site
# EfDiCMUN: community species effective richness of each site
# EfDiDOMI: dominant species effective richness of each site
# EfDiCOMM: common species effective richness of each site
# EfDiRARE: rare species effective richness of each site

rsb_CMUN <- cbind(com1[, 1:6], com1[,7:251]/com1[,252])
rsb_DOMI <- cbind(com1[, 1:6], com1[,7:251]/com1[,252])
rsb_COMM <- cbind(com1[, 1:6], com1[,7:251]/com1[,252])
rsb_RARE <- cbind(com1[, 1:6], com1[,7:251]/com1[,252])

for (i in 1:23) {
    # i = 1
    DOMI1 <- ifelse(spe_rsb[i, 7:251]>=5, 1, NA)
    COMM1 <- ifelse(spe_rsb[i, 7:251]>=1 & spe_rsb[i, 7:251]<5, 1, NA)
    RARE1 <- ifelse(spe_rsb[i, 7:251]<1, 1, NA)
    
    for (j in 7:251) {
        # j = 7
        rsb_DOMI[rsb_DOMI$Site == spe_rsb$Site[i], j] = rsb_DOMI[rsb_DOMI$Site == spe_rsb$Site[i], j] * DOMI1[j-6]
        rsb_COMM[rsb_COMM$Site == spe_rsb$Site[i], j] = rsb_COMM[rsb_COMM$Site == spe_rsb$Site[i], j] * COMM1[j-6]
        rsb_RARE[rsb_RARE$Site == spe_rsb$Site[i], j] = rsb_RARE[rsb_RARE$Site == spe_rsb$Site[i], j] * RARE1[j-6]
    }
}

richCMUN1 = apply(rsb_CMUN[, 7:251], 1, function(x)sum(ifelse(!is.na(x), 1, 0))); richCMUNI = ifelse(richCMUN1==0, NA, richCMUN1)
richDOMI1 = apply(rsb_DOMI[, 7:251], 1, function(x)sum(ifelse(!is.na(x), 1, 0))); richDOMII = ifelse(richDOMI1==0, NA, richDOMI1)
richCOMM1 = apply(rsb_COMM[, 7:251], 1, function(x)sum(ifelse(!is.na(x), 1, 0))); richCOMMI = ifelse(richCOMM1==0, NA, richCOMM1)
richRARE1 = apply(rsb_RARE[, 7:251], 1, function(x)sum(ifelse(!is.na(x), 1, 0))); richRAREI = ifelse(richRARE1==0, NA, richRARE1)

ShWiCMUN1 = apply(rsb_CMUN[, 7:251], 1, function(x)-sum(x*log(x), na.rm=T)); ShWiCMUNI = ifelse(ShWiCMUN1==0, NA, ShWiCMUN1)
ShWiDOMI1 = apply(rsb_DOMI[, 7:251], 1, function(x)-sum(x*log(x), na.rm=T)); ShWiDOMII = ifelse(ShWiDOMI1==0, NA, ShWiDOMI1)
ShWiCOMM1 = apply(rsb_COMM[, 7:251], 1, function(x)-sum(x*log(x), na.rm=T)); ShWiCOMMI = ifelse(ShWiCOMM1==0, NA, ShWiCOMM1)
ShWiRARE1 = apply(rsb_RARE[, 7:251], 1, function(x)-sum(x*log(x), na.rm=T)); ShWiRAREI = ifelse(ShWiRARE1==0, NA, ShWiRARE1)

richCMUN = apply(tapply(richCMUNI, list(rsb_CMUN$Site, rsb_CMUN$Year), mean, na.rm=T), 1, function(x)mean(x, na.rm=T))
richDOMI = apply(tapply(richDOMII, list(rsb_CMUN$Site, rsb_CMUN$Year), mean, na.rm=T), 1, function(x)mean(x, na.rm=T))
richCOMM = apply(tapply(richCOMMI, list(rsb_CMUN$Site, rsb_CMUN$Year), mean, na.rm=T), 1, function(x)mean(x, na.rm=T))
richRARE = apply(tapply(richRAREI, list(rsb_CMUN$Site, rsb_CMUN$Year), mean, na.rm=T), 1, function(x)mean(x, na.rm=T))

EfDiCMUN = exp(apply(tapply(ShWiCMUNI, list(rsb_CMUN$Site, rsb_CMUN$Year), mean, na.rm=T), 1, function(x)mean(x, na.rm=T)))
EfDiDOMI = exp(apply(tapply(ShWiDOMII, list(rsb_CMUN$Site, rsb_CMUN$Year), mean, na.rm=T), 1, function(x)mean(x, na.rm=T)))
EfDiCOMM = exp(apply(tapply(ShWiCOMMI, list(rsb_CMUN$Site, rsb_CMUN$Year), mean, na.rm=T), 1, function(x)mean(x, na.rm=T)))
EfDiRARE = exp(apply(tapply(ShWiRAREI, list(rsb_CMUN$Site, rsb_CMUN$Year), mean, na.rm=T), 1, function(x)mean(x, na.rm=T)))

################################################################################
# Following codes used to calculate:
# z_site: mean-variance scaling exponent of each site

fun_zscale = function(x, a, b) {
    a * x ^ b
}

z_par = matrix(NA, ncol=4, nrow=23)
for (i in 1:23) {
    # i = 1
    temp_spe_var1 = as.matrix(spe_var[i, 7:251])
    temp_spe_bm1 = as.matrix(spe_bm[i, 7:251])
    temp_spe_var = temp_spe_var1[which(!is.na(temp_spe_var1))]
    temp_spe_bm = temp_spe_bm1[which(!is.na(temp_spe_var1))]
    
    if (length(which(!is.na(temp_spe_var))) >= 5) {
        
        for (j in 1:5000) {
            temp_par = c(1, NA, 1.5, NA)
            
            nls_simu = nls(temp_spe_var~fun_zscale(temp_spe_bm, temp_par[1], b), start=list(b=temp_par[3]), na.action=na.omit)
            temp_par[3] = summary(nls_simu)[[10]][1, 1]
            temp_par[4] = summary(nls_simu)[[10]][1, 2]
            
            nls_simu = nls(temp_spe_var~fun_zscale(temp_spe_bm, a, temp_par[3]), start=list(a=temp_par[1]), na.action=na.omit)
            temp_par[1] = summary(nls_simu)[[10]][1, 1]
            temp_par[2] = summary(nls_simu)[[10]][1, 2]
        }
    } else {
        temp_par[1:4] = NA
    }
    
    z_par[i, ] = temp_par
}

a_site = z_par[, 1]
z_site = z_par[, 3]

################################################################################
# Following codes used to calculate:
# single_spe_tcv: species temporal CV of each species in each site
# single_spe_rsbm: relative species biomass of each species in each site

spe_rsbm1 = spe_rsb[, 7:251]
spe_sd1 = sqrt(spe_var[, 7:251])
spe_bm1 = NULL
for (i in 1:245) {
    spe_bm1 = cbind(spe_bm1, c(spe_rsbm1[, i] * spe_bm$cBM))
}
spe_tcv1 = NULL
for (i in 1:245) {
    spe_tcv1 = cbind(spe_tcv1, c(spe_sd1[, i] / spe_bm1[, i]))
    
}

single_spe_tcv = as.matrix(spe_tcv1); dim(single_spe_tcv) = c(245*23, 1)
single_spe_rsbm = as.matrix(spe_rsbm1); dim(single_spe_rsbm) = c(245*23, 1)

################################################################################
# Following codes used to calculate:
# sync_DD: synchronous species between dominant species
# sync_DC: synchronous species between dominant and common species
# sync_DR: synchronous species between dominant and rare species
# sync_CC: synchronous species between common species
# sync_CR: synchronous species between common and rare species
# sync_RR: synchronous species between rare species

sync_test = matrix(NA, ncol=7, nrow=23); sync_test
colnames(sync_test) = c("commu", "dom_dom", "dom_com", "dom_rare", "com_com", "com_rare", "rare_rare")
for (h in 1:23) {
    # h = 2
    temp1 = subset(spe_bm_year, Site==h)
    for (ti in 1:5) {
        # ti = 1
        if (all(is.na(temp1[ti, 7:251]))) {
            next
        } else {
            temp1[ti, 7:251] = ifelse(is.na(temp1[ti, 7:251]), 0, temp1[ti, 7:251])
        }
    }
    # names(temp1)
    # temp1[, 7]
    
    var_matrix = matrix(NA, ncol=245, nrow=245)
    for (i in 1:245) {
        for (j in 1:245) {
            var_matrix[i, j] = var(temp1[, 7:251][,i], temp1[, 7:251][,j], na.rm=T)
        }
    }
    
    dd_matrix = as.data.frame(var_matrix)
    for (ddi in 1:245) {
        dd_matrix[ddi, ] = dd_matrix[ddi, ] * spedom_bool[h, 7:251]
        dd_matrix[, ddi] = dd_matrix[, ddi] * t(spedom_bool[h, 7:251])
    }
    
    dc_matrix = as.data.frame(var_matrix)
    for (dci in 1:245) {
        dc_matrix[dci, ] = dc_matrix[dci, ] * spedom_bool[h, 7:251]
        dc_matrix[, dci] = dc_matrix[, dci] * t(specom_bool[h, 7:251])
    }
    
    dr_matrix = as.data.frame(var_matrix)
    for (dri in 1:245) {
        dr_matrix[dri, ] = dr_matrix[dri, ] * spedom_bool[h, 7:251]
        dr_matrix[, dri] = dr_matrix[, dri] * t(sperare_bool[h, 7:251])
    }
    
    cc_matrix = as.data.frame(var_matrix)
    for (cci in 1:245) {
        cc_matrix[cci, ] = cc_matrix[cci, ] * specom_bool[h, 7:251]
        cc_matrix[, cci] = cc_matrix[, cci] * t(specom_bool[h, 7:251])
    }
    
    cr_matrix = as.data.frame(var_matrix)
    for (cri in 1:245) {
        cr_matrix[cri, ] = cr_matrix[cri, ] * specom_bool[h, 7:251]
        cr_matrix[, cri] = cr_matrix[, cri] * t(sperare_bool[h, 7:251])
    }
    
    rr_matrix = as.data.frame(var_matrix)
    for (rri in 1:245) {
        rr_matrix[rri, ] = rr_matrix[rri, ] * sperare_bool[h, 7:251]
        rr_matrix[, rri] = rr_matrix[, rri] * t(sperare_bool[h, 7:251])
    }
    
    sync_test[h, 1] = sum(var_matrix, na.rm=T)/(sum(sqrt(diag(var_matrix)), na.rm=T)^2)
    sync_test[h, 2] = sum(dd_matrix, na.rm=T)/(sum(sqrt(diag(var_matrix)), na.rm=T)^2)
    sync_test[h, 3] = sum(dc_matrix, na.rm=T)/(sum(sqrt(diag(var_matrix)), na.rm=T)^2)
    sync_test[h, 4] = sum(dr_matrix, na.rm=T)/(sum(sqrt(diag(var_matrix)), na.rm=T)^2)
    sync_test[h, 5] = sum(cc_matrix, na.rm=T)/(sum(sqrt(diag(var_matrix)), na.rm=T)^2)
    sync_test[h, 6] = sum(cr_matrix, na.rm=T)/(sum(sqrt(diag(var_matrix)), na.rm=T)^2)
    sync_test[h, 7] = sum(rr_matrix, na.rm=T)/(sum(sqrt(diag(var_matrix)), na.rm=T)^2)
}

sync_DD = sync_test[, 2]
sync_DC = sync_test[, 3]
sync_DR = sync_test[, 4]
sync_CC = sync_test[, 5]
sync_CR = sync_test[, 6]
sync_RR = sync_test[, 7]

#########
# VegeType: Vegetation type of each site

VegeType = c("Meadow steppe", "Meadow steppe", "Typical steppe", "Typical steppe",
             "Typical steppe", "Typical steppe", "Typical steppe", "Typical steppe",
             "Typical steppe", "Meadow steppe", "Typical steppe", "Typical steppe",
             "Typical steppe", "Typical steppe", "Typical steppe", "Typical steppe",
             "Typical steppe", "Desert steppe", "Desert steppe", "Desert steppe",
             "Desert steppe", "Typical steppe", "Desert steppe")

########
# Output variables

data_table = cbind(spe_rsb[,2:4], VegeType, pre_site, pre_cv, bm_site, 
                   richCMUN, richDOMI, richCOMM, richRARE, 
                   EfDiCMUN, EfDiDOMI, EfDiCOMM, EfDiRARE, 
                   cbm_cv, w_specv_site, sync_site, z_site, 
                   w_spedomcv_site, w_specomcv_site, w_sperarecv_site, dom_spesd_site, dom_spebm_site,
                   sync_DD, sync_DC, sync_DR, sync_CC, sync_CR, sync_RR)

colnames(data_table) = c("Site", "Latitude", "Longitude", "VegeType", "MGP", "CV_MGP", "CMU_BM",
                         "CMU_SR", "D_SR", "C_SR", "R_SR",
                         "CMU_ER", "D_ER", "C_ER", "R_ER",
                         "CMU_CV", "CMU_SCV", "CMU_SYNC", "Taylor_Z",
                         "D_SCV", "C_SCV", "R_SCV", "DOM_SD", "DOM_BM",
                         "DD_SYNC", "DC_SYNC", "DR_SYNC", "CC_SYNC", "CR_SYNC", "RR_SYNC")

# write.csv(data_table, "Wang_et_al._2020_Dataset.csv", row.names=F)
