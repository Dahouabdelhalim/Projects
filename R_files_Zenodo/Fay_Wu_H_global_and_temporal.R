#load the required library 
library(ggplot2)
library(PopGenome)

#read in the vcf files for the two loci (divergent and neutral); these were obtained using the "make_input_files_for_Fay_Wu.chr7_locus.sh" script
divergent_locus<-readVCF(filename = "outgroup_and_sagrei.chr7.divergent_locus.vcf.gz",frompos = 70000000,topos = 88000000, tid = "7",numcols = 500)
neutral_locus<-readVCF(filename = "outgroup_and_sagrei.chr7.neutral_locus.vcf.gz",frompos = 21000000,topos = 59114159, tid = "7",numcols = 500)

#set the populations (the 2018 and 2003 timepoints for invasive populations are treated as separate populations) and the outgroup sequence for both datasets
divergent_locus<-set.populations(divergent_locus,list(c("ALA_D1","ALA_D2","ALA_D3","ALA_D4","ALA_D5","ALA_D6","ALA_J1","ALA_J2","ALA_J4","ALA_J5","ALA_J6","ALA_Z1","ALA_Z2","ALA_Z4","ALA_Z5","ALA_Z6","ALA_Z7"),
                                                      c("GVL_1301","GVL_1305","GVL_1309","GVL_1312","GVL_1313","GVL_1314","GVL_1315"),
                                                c("BRE_D2","BRE_D3","BRE_D4","BRE_D5","BRE_D6","BRE_D7R","BRE_J1","BRE_J2","BRE_J3","BRE_Z1","BRE_Z10","BRE_Z2","BRE_Z3","BRE_Z4","BRE_Z5","BRE_Z6","BRE_Z7","BRE_Z8","BRE_Z9"),
                                                c("BRO_D1","BRO_D2","BRO_D3","BRO_D4","BRO_D5","BRO_D6","BRO_D7","BRO_J1","BRO_J2","BRO_Z1","BRO_Z10","BRO_Z11","BRO_Z2","BRO_Z3p6","BRO_Z4","BRO_Z5","BRO_Z6","BRO_Z7","BRO_Z8","BRO_Z9"),
                                                c("CIT_D1","CIT_J2","CIT_J3","CIT_J4","CIT_J5","CIT_J6","CIT_S1","CIT_S10","CIT_S2","CIT_S3","CIT_S4","CIT_S5","CIT_S6","CIT_S7","CIT_S8","CIT_S9"),
                                                c("CMB_D1","CMB_D2","CMB_D3","CMB_D4","CMB_D5","CMB_D7","CMB_D8","CMB_J1","CMB_J2","CMB_J3","CMB_J4","CMB_J5","CMB_J6","CMB_Z1","CMB_Z2R","CMB_Z3","CMB_Z4","CMB_Z6"),
                                                c("COL_D1","COL_D2","COL_J1","COL_J2","COL_J3","COL_J4","COL_J5","COL_S1","COL_S2","COL_S3","COL_S4","COL_S5","COL_S6","COL_S7","COL_S8"),
                                                c("DUV_D1","DUV_D2","DUV_D3","DUV_D4","DUV_D5","DUV_D6","DUV_D7","DUV_J1","DUV_J2","DUV_Z1","DUV_Z10","DUV_Z11","DUV_Z3","DUV_Z4","DUV_Z5","DUV_Z6","DUV_Z8","DUV_Z9"),
                                                c("FLA_D1","FLA_D2","FLA_D3","FLA_D4","FLA_J1","FLA_J2","FLA_Z1","FLA_Z10","FLA_Z10_2","FLA_Z11","FLA_Z12","FLA_Z13","FLA_Z2","FLA_Z3","FLA_Z4","FLA_Z6","FLA_Z7","FLA_Z8","FLA_Z9"),
                                                c("GLA_D2","GLA_D3","GLA_D4","GLA_D5","GLA_J1","GLA_J2","GLA_J3","GLA_J4","GLA_Z1","GLA_Z10","GLA_Z11","GLA_Z2","GLA_Z3","GLA_Z4","GLA_Z5","GLA_Z6","GLA_Z7","GLA_Z8","GLA_Z9"),
                                                c("HAM_D1","HAM_D10","HAM_D2","HAM_D3","HAM_D4","HAM_D5","HAM_D7","HAM_D8","HAM_D9","HAM_J1","HAM_J2","HAM_J3","HAM_J4","HAM_J5","HAM_J6","HAM_R1","HAM_R10","HAM_R2","HAM_R3","HAM_R4","HAM_R5","HAM_R6","HAM_R7","HAM_R8","HAM_R9"),
                                                c("HEN_D1","HEN_D2","HEN_D3","HEN_D4","HEN_D5","HEN_D6","HEN_D7","HEN_J1","HEN_J2","HEN_Z1","HEN_Z10","HEN_Z11","HEN_Z2","HEN_Z3","HEN_Z5","HEN_Z6","HEN_Z7","HEN_Z8"),
                                                c("HIG_D1","HIG_D2","HIG_J1","HIG_J2","HIG_J4","HIG_Z1","HIG_Z10","HIG_Z11","HIG_Z12","HIG_Z13","HIG_Z2","HIG_Z3","HIG_Z4","HIG_Z5","HIG_Z6","HIG_Z7","HIG_Z8","HIG_Z9"),
                                                c("LAK_D1","LAK_D2","LAK_D3","LAK_D4","LAK_D5","LAK_D6","LAK_J1","LAK_J2","LAK_J3","LAK_J4","LAK_J5","LAK_J6","LAK_Z1","LAK_Z2","LAK_Z3","LAK_Z4","LAK_Z5","LAK_Z6","LAK_Z7","LAK_Z8"),
                                                c("LEE_D1","LEE_D2","LEE_D3","LEE_D4","LEE_D5","LEE_J1","LEE_J2","LEE_R2","LEE_R3","LEE_R4","LEE_R5","LEE_R8","LEE_R9","LEE_Runknown","LEE_S1","LEE_S2","LEE_S3","LEE_S4","LEE_S5","LEE_S6","LEE_S7","LEE_S8","LEE_S9"),
                                                c("LEV_D1","LEV_D2","LEV_D3","LEV_D4","LEV_D5","LEV_D6","LEV_D7","LEV_D8","LEV_J1","LEV_J2","LEV_J3","LEV_J4","LEV_J5","LEV_J6","LEV_J7","LEV_J8","LEV_R2","LEV_R3","LEV_R4","LEV_R5","LEV_R6","LEV_R7","LEV_R8","LEV_R9"),
                                                c("LOW_1","LOW_10","LOW_11","LOW_12","LOW_13","LOW_14","LOW_15","LOW_16","LOW_17","LOW_18","LOW_19","LOW_2","LOW_20","LOW_3","LOW_4","LOW_5","LOW_6","LOW_7","LOW_8","LOW_9"),
                                                c("MAN_D1","MAN_D2","MAN_D3","MAN_D4","MAN_D5","MAN_D6","MAN_J1","MAN_J10","MAN_J2","MAN_J3","MAN_J4","MAN_J5","MAN_J6","MAN_J7","MAN_J8","MAN_J9","MAN_R1","MAN_R3","MAN_R4","MAN_R6","MAN_R7","MAN_R8","MAN_X1","MAN_X2","MAN_X3","MAN_X4"),
                                                c("MAR_D1","MAR_D2","MAR_D3","MAR_D4","MAR_D5","MAR_J1","MAR_J2","MAR_J3","MAR_J5","MAR_J6","MAR_J7","MAR_Z1","MAR_Z2","MAR_Z3","MAR_Z4","MAR_Z5","MAR_Z6","MAR_Z7","MAR_Z8"),
                                                c("MIA_D1","MIA_D2","MIA_D3","MIA_D4","MIA_D5","MIA_D6","MIA_D7","MIA_J1","MIA_J2","MIA_Z1","MIA_Z10","MIA_Z2","MIA_Z3","MIA_Z4","MIA_Z5","MIA_Z6","MIA_Z7","MIA_Z8","MIA_Z9"),
                                                c("CGA_406","CGA_408","CGA_409","CGA_410","CGA_453","CGA_455","CGA_456","CGAq_407"),
                                                c("MON_D2","MON_D3","MON_D4","MON_D5","MON_D6","MON_J1","MON_J2","MON_J3","MON_Z1","MON_Z10","MON_Z11","MON_Z2","MON_Z3","MON_Z5","MON_Z6","MON_Z7","MON_Z8","MON_Z9"),
                                                c("LMA_376","LMA_377","LMA_378","LMA_379","LMA_426","LMA_427","LMA_428","LMA_429","LMA_430","LMA_437"),
                                                c("ORA_D1","ORA_D2","ORA_D3","ORA_D4","ORA_D5","ORA_J1","ORA_J2","ORA_J4","ORA_J5","ORA_Z10","ORA_Z2","ORA_Z3","ORA_Z4","ORA_Z6","ORA_Z7","ORA_Z8","ORA_Z9"),
                                                c("ORL_493","ORL_494","ORL_496","ORL_497","ORL_500","ORL_506","ORL_511"),
                                                c("PAL_D1","PAL_D10","PAL_D2","PAL_D3","PAL_D4","PAL_D5","PAL_D6","PAL_D7","PAL_D8","PAL_D9","PAL_J1","PAL_J2","PAL_J3","PAL_Z1","PAL_Z2","PAL_Z3","PAL_Z4","PAL_Z5","PAL_Z6","PAL_Z7"),
                                                c("POL_D1","POL_D2","POL_D3","POL_D4","POL_J1","POL_J2","POL_J3","POL_J4","POL_J5","POL_J6","POL_J7","POL_Z1","POL_Z3","POL_Z4","POL_Z5","POL_Z6","POL_Z7","POL_Z8","POL_Z9"),
                                                c("SAR_D1","SAR_D2","SAR_D3","SAR_J1","SAR_J2","SAR_J3","SAR_J4","SAR_R1","SAR_R2","SAR_R3","SAR_R4","SAR_R6","SAR_R7","SAR_R8","SAR_S1","SAR_S2","SAR_S3","SAR_S4","SAR_S5","SAR_S6","SAR_S7","SAR_S8","SAR_S9"),
                                                c("STJ_D1","STJ_D2","STJ_J1","STJ_J2","STJ_Z1","STJ_Z2","STJ_Z4","STJ_Z5","STJ_Z6","STJ_Z7","STJ_Z8","STJ_Z9"),
                                                c("STL_D1","STL_D2","STL_D3","STL_D4","STL_D5","STL_J1","STL_J2","STL_J3","STL_J4","STL_Z1","STL_Z2","STL_Z3","STL_Z4","STL_Z5","STL_Z6","STL_Z7"),
                                                c("STP_D1","STP_D4R","STP_D5","STP_D6","STP_J1","STP_J3","STP_S1R","STP_S2","STP_S3","STP_S4","STP_S5","STP_S6","STP_S7"),
                                                c("STP_282","STP_283","STP_285","STP_286","STP_323","STP_324","STP_325","STP_326","STP_327","STP_329","STP_330","STP_331","STP_332","STP_333","STP_334","STP_335","STP_336","STP_337"),
                                                c("TAM_D1","TAM_D2","TAM_J1","TAM_J2","TAM_J3","TAM_S1","TAM_S10","TAM_S11","TAM_S2","TAM_S3","TAM_S4","TAM_S5","TAM_S6","TAM_S7","TAM_S8","TAM_S9","TAM_X1","TAM_X2","TAM_X3","TAM_X4"),
                                                c("TAM_201","TAM_202","TAM_203","TAM_204","TAM_205","TAM_206","TAM_207","TAM_208","TAM_230","TAM_231","TAM_232","TAM_233","TAM_234","TAM_236","TAM_237","TAM_238","TAM_257","TAM_258"),
                                                c("TIF_1","TIF_10","TIF_11","TIF_12","TIF_13","TIF_14","TIF_15","TIF_16","TIF_17","TIF_18","TIF_2","TIF_3","TIF_4","TIF_5","TIF_6","TIF_7","TIF_8","TIF_9"),
                                                c("VOL_D1","VOL_D2","VOL_D4","VOL_D5","VOL_J1","VOL_Z1","VOL_Z4","VOL_Z5","VOL_Z6","VOL_Z8","VOL_Z9"),
                                                c("CAB_2918","CAB_2919","CAB_2920","CAB_2922","CAB_2925","CAB_2927","CAB_2930","CAB_2931"),
                                                c("ESM_2881","ESM_2882","ESM_2883","ESM_2884","ESM_2886","ESM_2887","ESM_2890","ESM_2893","ESM_2894","ESM_2896","ESM_2909","ESM_2910"),
                                                c("GUA_3050","GUA_3052","GUA_3053","GUA_3054","GUA_3055","GUA_3056","GUA_3057","GUA_3058","GUA_3062"),
                                                c("MAR_2971","MAR_2972","MAR_2973","MAR_2975","MAR_2977","MAR_2978","MAR_2979","MAR_2980","MAR_2981","MAR_2990","MAR_2993","MAR_2994","MAR_2995","MAR_2997","MAR_2998","MAR_2999","MAR_3000"),
                                                c("JIC_2932","JIC_2933","JIC_2934","JIC_2935","JIC_2936","JIC_2937","JIC_2938","JIC_2940","JIC_2941","JIC_2942","JIC_2943","JIC_2944","JIC_2945","JIC_2946","JIC_2947","JIC_2958","JIC_2959","JIC_2960","JIC_2962"),
                                                c("POR_02","POR_03","POR_04","POR_05","POR_06","POR_07","POR_09","POR_10","POR_11","POR_12","POR_13","POR_14","POR_15","POR_16","POR_17","POR_18","POR_19","POR_20","POR_21","POR_22"),
                                                c("sagrei_157","sagrei_159","sagrei_160","sagrei_161","sagrei_162","sagrei_163","sagrei_172","sagrei_173","sagrei_174","sagrei_175","sagrei_176","sagrei_177","sagrei_178","sagrei_179"),
                                                c("SOR_3004","SOR_3005","SOR_3006","SOR_3007","SOR_3008","SOR_3009","SOR_3010","SOR_3016","SOR_3017","SOR_3019","SOR_3021","SOR_3022"),
                                                c("LHB_2836","LHB_2837","LHB_2838","LHB_2839","LHB_2840","LHB_2841","LHB_2842","LHB_2843","LHB_2844","LHB_2845","LHB_2846","LHB_2849","LHB_2853","LHB_2854","LHB_2855"),
                                                c("LHA_2802","LHA_2813","LHA_2814","LHA_2815","LHA_2816","LHA_2817","LHA_2818","LHA_2820")))
divergent_locus <- set.outgroup(divergent_locus,c("SRR_8983021"))


neutral_locus<-set.populations(neutral_locus,list(c("ALA_D1","ALA_D2","ALA_D3","ALA_D4","ALA_D5","ALA_D6","ALA_J1","ALA_J2","ALA_J4","ALA_J5","ALA_J6","ALA_Z1","ALA_Z2","ALA_Z4","ALA_Z5","ALA_Z6","ALA_Z7"),
                                                      c("GVL_1301","GVL_1305","GVL_1309","GVL_1312","GVL_1313","GVL_1314","GVL_1315"),
                                                      c("BRE_D2","BRE_D3","BRE_D4","BRE_D5","BRE_D6","BRE_D7R","BRE_J1","BRE_J2","BRE_J3","BRE_Z1","BRE_Z10","BRE_Z2","BRE_Z3","BRE_Z4","BRE_Z5","BRE_Z6","BRE_Z7","BRE_Z8","BRE_Z9"),
                                                      c("BRO_D1","BRO_D2","BRO_D3","BRO_D4","BRO_D5","BRO_D6","BRO_D7","BRO_J1","BRO_J2","BRO_Z1","BRO_Z10","BRO_Z11","BRO_Z2","BRO_Z3p6","BRO_Z4","BRO_Z5","BRO_Z6","BRO_Z7","BRO_Z8","BRO_Z9"),
                                                      c("CIT_D1","CIT_J2","CIT_J3","CIT_J4","CIT_J5","CIT_J6","CIT_S1","CIT_S10","CIT_S2","CIT_S3","CIT_S4","CIT_S5","CIT_S6","CIT_S7","CIT_S8","CIT_S9"),
                                                      c("CMB_D1","CMB_D2","CMB_D3","CMB_D4","CMB_D5","CMB_D7","CMB_D8","CMB_J1","CMB_J2","CMB_J3","CMB_J4","CMB_J5","CMB_J6","CMB_Z1","CMB_Z2R","CMB_Z3","CMB_Z4","CMB_Z6"),
                                                      c("COL_D1","COL_D2","COL_J1","COL_J2","COL_J3","COL_J4","COL_J5","COL_S1","COL_S2","COL_S3","COL_S4","COL_S5","COL_S6","COL_S7","COL_S8"),
                                                      c("DUV_D1","DUV_D2","DUV_D3","DUV_D4","DUV_D5","DUV_D6","DUV_D7","DUV_J1","DUV_J2","DUV_Z1","DUV_Z10","DUV_Z11","DUV_Z3","DUV_Z4","DUV_Z5","DUV_Z6","DUV_Z8","DUV_Z9"),
                                                      c("FLA_D1","FLA_D2","FLA_D3","FLA_D4","FLA_J1","FLA_J2","FLA_Z1","FLA_Z10","FLA_Z10_2","FLA_Z11","FLA_Z12","FLA_Z13","FLA_Z2","FLA_Z3","FLA_Z4","FLA_Z6","FLA_Z7","FLA_Z8","FLA_Z9"),
                                                      c("GLA_D2","GLA_D3","GLA_D4","GLA_D5","GLA_J1","GLA_J2","GLA_J3","GLA_J4","GLA_Z1","GLA_Z10","GLA_Z11","GLA_Z2","GLA_Z3","GLA_Z4","GLA_Z5","GLA_Z6","GLA_Z7","GLA_Z8","GLA_Z9"),
                                                      c("HAM_D1","HAM_D10","HAM_D2","HAM_D3","HAM_D4","HAM_D5","HAM_D7","HAM_D8","HAM_D9","HAM_J1","HAM_J2","HAM_J3","HAM_J4","HAM_J5","HAM_J6","HAM_R1","HAM_R10","HAM_R2","HAM_R3","HAM_R4","HAM_R5","HAM_R6","HAM_R7","HAM_R8","HAM_R9"),
                                                      c("HEN_D1","HEN_D2","HEN_D3","HEN_D4","HEN_D5","HEN_D6","HEN_D7","HEN_J1","HEN_J2","HEN_Z1","HEN_Z10","HEN_Z11","HEN_Z2","HEN_Z3","HEN_Z5","HEN_Z6","HEN_Z7","HEN_Z8"),
                                                      c("HIG_D1","HIG_D2","HIG_J1","HIG_J2","HIG_J4","HIG_Z1","HIG_Z10","HIG_Z11","HIG_Z12","HIG_Z13","HIG_Z2","HIG_Z3","HIG_Z4","HIG_Z5","HIG_Z6","HIG_Z7","HIG_Z8","HIG_Z9"),
                                                      c("LAK_D1","LAK_D2","LAK_D3","LAK_D4","LAK_D5","LAK_D6","LAK_J1","LAK_J2","LAK_J3","LAK_J4","LAK_J5","LAK_J6","LAK_Z1","LAK_Z2","LAK_Z3","LAK_Z4","LAK_Z5","LAK_Z6","LAK_Z7","LAK_Z8"),
                                                      c("LEE_D1","LEE_D2","LEE_D3","LEE_D4","LEE_D5","LEE_J1","LEE_J2","LEE_R2","LEE_R3","LEE_R4","LEE_R5","LEE_R8","LEE_R9","LEE_Runknown","LEE_S1","LEE_S2","LEE_S3","LEE_S4","LEE_S5","LEE_S6","LEE_S7","LEE_S8","LEE_S9"),
                                                      c("LEV_D1","LEV_D2","LEV_D3","LEV_D4","LEV_D5","LEV_D6","LEV_D7","LEV_D8","LEV_J1","LEV_J2","LEV_J3","LEV_J4","LEV_J5","LEV_J6","LEV_J7","LEV_J8","LEV_R2","LEV_R3","LEV_R4","LEV_R5","LEV_R6","LEV_R7","LEV_R8","LEV_R9"),
                                                      c("LOW_1","LOW_10","LOW_11","LOW_12","LOW_13","LOW_14","LOW_15","LOW_16","LOW_17","LOW_18","LOW_19","LOW_2","LOW_20","LOW_3","LOW_4","LOW_5","LOW_6","LOW_7","LOW_8","LOW_9"),
                                                      c("MAN_D1","MAN_D2","MAN_D3","MAN_D4","MAN_D5","MAN_D6","MAN_J1","MAN_J10","MAN_J2","MAN_J3","MAN_J4","MAN_J5","MAN_J6","MAN_J7","MAN_J8","MAN_J9","MAN_R1","MAN_R3","MAN_R4","MAN_R6","MAN_R7","MAN_R8","MAN_X1","MAN_X2","MAN_X3","MAN_X4"),
                                                      c("MAR_D1","MAR_D2","MAR_D3","MAR_D4","MAR_D5","MAR_J1","MAR_J2","MAR_J3","MAR_J5","MAR_J6","MAR_J7","MAR_Z1","MAR_Z2","MAR_Z3","MAR_Z4","MAR_Z5","MAR_Z6","MAR_Z7","MAR_Z8"),
                                                      c("MIA_D1","MIA_D2","MIA_D3","MIA_D4","MIA_D5","MIA_D6","MIA_D7","MIA_J1","MIA_J2","MIA_Z1","MIA_Z10","MIA_Z2","MIA_Z3","MIA_Z4","MIA_Z5","MIA_Z6","MIA_Z7","MIA_Z8","MIA_Z9"),
                                                      c("CGA_406","CGA_408","CGA_409","CGA_410","CGA_453","CGA_455","CGA_456","CGAq_407"),
                                                      c("MON_D2","MON_D3","MON_D4","MON_D5","MON_D6","MON_J1","MON_J2","MON_J3","MON_Z1","MON_Z10","MON_Z11","MON_Z2","MON_Z3","MON_Z5","MON_Z6","MON_Z7","MON_Z8","MON_Z9"),
                                                      c("LMA_376","LMA_377","LMA_378","LMA_379","LMA_426","LMA_427","LMA_428","LMA_429","LMA_430","LMA_437"),
                                                      c("ORA_D1","ORA_D2","ORA_D3","ORA_D4","ORA_D5","ORA_J1","ORA_J2","ORA_J4","ORA_J5","ORA_Z10","ORA_Z2","ORA_Z3","ORA_Z4","ORA_Z6","ORA_Z7","ORA_Z8","ORA_Z9"),
                                                      c("ORL_493","ORL_494","ORL_496","ORL_497","ORL_500","ORL_506","ORL_511"),
                                                      c("PAL_D1","PAL_D10","PAL_D2","PAL_D3","PAL_D4","PAL_D5","PAL_D6","PAL_D7","PAL_D8","PAL_D9","PAL_J1","PAL_J2","PAL_J3","PAL_Z1","PAL_Z2","PAL_Z3","PAL_Z4","PAL_Z5","PAL_Z6","PAL_Z7"),
                                                      c("POL_D1","POL_D2","POL_D3","POL_D4","POL_J1","POL_J2","POL_J3","POL_J4","POL_J5","POL_J6","POL_J7","POL_Z1","POL_Z3","POL_Z4","POL_Z5","POL_Z6","POL_Z7","POL_Z8","POL_Z9"),
                                                      c("SAR_D1","SAR_D2","SAR_D3","SAR_J1","SAR_J2","SAR_J3","SAR_J4","SAR_R1","SAR_R2","SAR_R3","SAR_R4","SAR_R6","SAR_R7","SAR_R8","SAR_S1","SAR_S2","SAR_S3","SAR_S4","SAR_S5","SAR_S6","SAR_S7","SAR_S8","SAR_S9"),
                                                      c("STJ_D1","STJ_D2","STJ_J1","STJ_J2","STJ_Z1","STJ_Z2","STJ_Z4","STJ_Z5","STJ_Z6","STJ_Z7","STJ_Z8","STJ_Z9"),
                                                      c("STL_D1","STL_D2","STL_D3","STL_D4","STL_D5","STL_J1","STL_J2","STL_J3","STL_J4","STL_Z1","STL_Z2","STL_Z3","STL_Z4","STL_Z5","STL_Z6","STL_Z7"),
                                                      c("STP_D1","STP_D4R","STP_D5","STP_D6","STP_J1","STP_J3","STP_S1R","STP_S2","STP_S3","STP_S4","STP_S5","STP_S6","STP_S7"),
                                                      c("STP_282","STP_283","STP_285","STP_286","STP_323","STP_324","STP_325","STP_326","STP_327","STP_329","STP_330","STP_331","STP_332","STP_333","STP_334","STP_335","STP_336","STP_337"),
                                                      c("TAM_D1","TAM_D2","TAM_J1","TAM_J2","TAM_J3","TAM_S1","TAM_S10","TAM_S11","TAM_S2","TAM_S3","TAM_S4","TAM_S5","TAM_S6","TAM_S7","TAM_S8","TAM_S9","TAM_X1","TAM_X2","TAM_X3","TAM_X4"),
                                                      c("TAM_201","TAM_202","TAM_203","TAM_204","TAM_205","TAM_206","TAM_207","TAM_208","TAM_230","TAM_231","TAM_232","TAM_233","TAM_234","TAM_236","TAM_237","TAM_238","TAM_257","TAM_258"),
                                                      c("TIF_1","TIF_10","TIF_11","TIF_12","TIF_13","TIF_14","TIF_15","TIF_16","TIF_17","TIF_18","TIF_2","TIF_3","TIF_4","TIF_5","TIF_6","TIF_7","TIF_8","TIF_9"),
                                                      c("VOL_D1","VOL_D2","VOL_D4","VOL_D5","VOL_J1","VOL_Z1","VOL_Z4","VOL_Z5","VOL_Z6","VOL_Z8","VOL_Z9"),
                                                      c("CAB_2918","CAB_2919","CAB_2920","CAB_2922","CAB_2925","CAB_2927","CAB_2930","CAB_2931"),
                                                      c("ESM_2881","ESM_2882","ESM_2883","ESM_2884","ESM_2886","ESM_2887","ESM_2890","ESM_2893","ESM_2894","ESM_2896","ESM_2909","ESM_2910"),
                                                      c("GUA_3050","GUA_3052","GUA_3053","GUA_3054","GUA_3055","GUA_3056","GUA_3057","GUA_3058","GUA_3062"),
                                                      c("MAR_2971","MAR_2972","MAR_2973","MAR_2975","MAR_2977","MAR_2978","MAR_2979","MAR_2980","MAR_2981","MAR_2990","MAR_2993","MAR_2994","MAR_2995","MAR_2997","MAR_2998","MAR_2999","MAR_3000"),
                                                      c("JIC_2932","JIC_2933","JIC_2934","JIC_2935","JIC_2936","JIC_2937","JIC_2938","JIC_2940","JIC_2941","JIC_2942","JIC_2943","JIC_2944","JIC_2945","JIC_2946","JIC_2947","JIC_2958","JIC_2959","JIC_2960","JIC_2962"),
                                                      c("POR_02","POR_03","POR_04","POR_05","POR_06","POR_07","POR_09","POR_10","POR_11","POR_12","POR_13","POR_14","POR_15","POR_16","POR_17","POR_18","POR_19","POR_20","POR_21","POR_22"),
                                                      c("sagrei_157","sagrei_159","sagrei_160","sagrei_161","sagrei_162","sagrei_163","sagrei_172","sagrei_173","sagrei_174","sagrei_175","sagrei_176","sagrei_177","sagrei_178","sagrei_179"),
                                                      c("SOR_3004","SOR_3005","SOR_3006","SOR_3007","SOR_3008","SOR_3009","SOR_3010","SOR_3016","SOR_3017","SOR_3019","SOR_3021","SOR_3022"),
                                                      c("LHB_2836","LHB_2837","LHB_2838","LHB_2839","LHB_2840","LHB_2841","LHB_2842","LHB_2843","LHB_2844","LHB_2845","LHB_2846","LHB_2849","LHB_2853","LHB_2854","LHB_2855"),
                                                      c("LHA_2802","LHA_2813","LHA_2814","LHA_2815","LHA_2816","LHA_2817","LHA_2818","LHA_2820")))
neutral_locus <- set.outgroup(neutral_locus,c("SRR_8983021"))

#population names, used for the output H index
populations<-c("ALA","ALA2003","BRE","BRO","CIT","CMB","COL","DUV","FLA","GLA","HAM","HEN","HIG","LAK","LEE","LEV","LOW","MAN","MARi","MIA","MIA2003","MON","MON2003","ORA","ORA2003","PAL","POL","SAR","STJ","STL","STP","STP2003","TAM","TAM2003","TIF","VOL","CAB","ESM","GUA","MAR","JIC","POR","BAH","SOR","LHB","LHA")

#compute neutrality statistics and get the H index values
divergent_locus  <- neutrality.stats(divergent_locus)
get.neutrality(divergent_locus)
H_index_divergent<-as.data.frame(divergent_locus@Fay.Wu.H)
names(H_index_divergent) <- populations
write.table(H_index_divergent, file="H_index_divergent_locus_with_2003_samples.csv",sep = '\\t', row.names = FALSE)

neutral_locus  <- neutrality.stats(neutral_locus)
get.neutrality(neutral_locus)
H_index_neutral<-as.data.frame(neutral_locus@Fay.Wu.H)
names(H_index_neutral) <- populations
write.table(H_index_neutral, file="H_index_neutral_locus_with_2003_samples.csv",sep = '\\t', row.names = FALSE)


###the input file contains H results obtained above (H_index_divergent_locus_with_2003_samples.csv and H_index_neutral_locus_with_2003_samples.csv), 
#to which we added information on locus category (divergent or neutral), the range of each population (Invasive or Native), year (2003 or other), and whether the H estimates will be used in the temporal analysis (6 invasive populations, 2003 vs. 2018 invasive populations).  
H_results<-read.csv(file="H_index_results.csv", header = TRUE)
H_results_global<- H_results[H_results$Year != '2003',]
H_results_temporal<- H_results[H_results$Analysis == 'temporal_analysis',]

#global analysis first (native populations and 2018 invasive populations)
#make data subsets for the global analysis
native_H_results <- H_results_global[H_results_global$Range == 'Native',]
native_divergent_H_results <- native_H_results[native_H_results$Locus_category == 'divergent',]
native_neutral_H_results <- native_H_results[native_H_results$Locus_category == 'neutral',]

invasive_H_results <- H_results_global[H_results_global$Range == 'Invasive',]
invasive_divergent_H_results <- invasive_H_results[invasive_H_results$Locus_category == 'divergent',]
invasive_neutral_H_results <- invasive_H_results[invasive_H_results$Locus_category == 'neutral',]

#test whether divergent locus H scores are smaller than neutral locus H scores, for each range
wilcox.test(native_divergent_H_results$H, native_neutral_H_results$H,alternative = "less", paired = TRUE)
wilcox.test(invasive_divergent_H_results$H, invasive_neutral_H_results$H,alternative = "less", paired = TRUE)

#test whether native range H scores are smaller than invasive range H scores, for each locus category
wilcox.test(native_divergent_H_results$H, invasive_divergent_H_results$H,alternative = "less")
wilcox.test(native_neutral_H_results$H, invasive_neutral_H_results$H,alternative = "less")

#adjust P values for multiple comparisons using Bonferroni
pvalues<-c("0.004883","0.8973","3.173e-07","0.8423")
p.adjust(pvalues, method="bonferroni")


#temporal analysis (6 invasive populations, 2003 vs. 2018 invasive populations)
#make data subsets for the temporal analysis
H_results_temporal_2003 <- H_results_temporal[H_results_temporal$Year == '2003',]
H_results_temporal_2018 <- H_results_temporal[H_results_temporal$Year != '2003',]

H_results_temporal_2003_divergent <- H_results_temporal_2003[H_results_temporal_2003$Locus_category == 'divergent',]
H_results_temporal_2003_neutral <- H_results_temporal_2003[H_results_temporal_2003$Locus_category == 'neutral',]

H_results_temporal_2018_divergent <- H_results_temporal_2018[H_results_temporal_2018$Locus_category == 'divergent',]
H_results_temporal_2018_neutral <- H_results_temporal_2018[H_results_temporal_2018$Locus_category == 'neutral',]

#test whether divergent locus H scores are smaller than neutral locus H scores, for each timepoint
wilcox.test(H_results_temporal_2003_divergent$H, H_results_temporal_2003_neutral$H,alternative = "less", paired = TRUE)
wilcox.test(H_results_temporal_2018_divergent$H, H_results_temporal_2018_neutral$H,alternative = "less", paired = TRUE)

#test whether divergent locus H scores are different between the two timepoints
wilcox.test(H_results_temporal_2003_divergent$H, H_results_temporal_2018_divergent$H, paired = TRUE)

#adjust P values for multiple comparisons using Bonferroni
pvalues<-c("0.8906","0.2188","0.1563")
p.adjust(pvalues, method="bonferroni")


#set the order of locus categories (neutral then divergent), and plot results (as per Fig. 2D)
native_H_results$Locus_category = factor(native_H_results$Locus_category, levels=c("neutral", "divergent"))
invasive_H_results$Locus_category = factor(invasive_H_results$Locus_category, levels=c("neutral", "divergent"))

#native range H plot
ggplot(data = native_H_results, aes(x=Locus_category, y=H))+
  geom_line(aes(group = Population_ID))+
  geom_point(size=5, col="red")+
  labs(x = "Locus", y = "Fay & Wu's H") + 
  ylim(-4,1.5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        axis.text=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

#invasive range H plot
ggplot(data = invasive_H_results, aes(x=Locus_category, y=H))+
  geom_line(aes(group = Population_ID))+
  geom_point(size=5, col="black")+
  labs(x = "Locus", y = "Fay & Wu's H") + 
  ylim(-4,1.5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        axis.text=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


#2003 temporal analysis H plot
H_results_temporal_2003$Locus_category = factor(H_results_temporal_2003$Locus_category, levels=c("neutral", "divergent"))
H_results_temporal_2018$Locus_category = factor(H_results_temporal_2018$Locus_category, levels=c("neutral", "divergent"))

ggplot(data = H_results_temporal_2003, aes(x=Locus_category, y=H))+
  geom_line(aes(group = Population_ID))+
  geom_point(size=5)+
  labs(x = "Locus", y = "Fay & Wu's H") + 
  ylim(-4,1.5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        axis.text=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

#2018 temporal analysis H plot
ggplot(data = H_results_temporal_2018, aes(x=Locus_category, y=H))+
  geom_line(aes(group = Population_ID))+
  geom_point(size=5)+
  labs(x = "Locus", y = "Fay & Wu's H") + 
  ylim(-4,1.5)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        axis.text=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())




