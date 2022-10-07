Signature optimization with a genetic algorithm
================

This Markdown illustrates the methodology for optimization the COVID-19
signature using a genetic algorithm. The optimization parameters
(population, max iterations, and grid size) are set to minimize the
computational time when generating the Markdown file. Note that, when
using more extensive searches, the computational time increases
substantially and may require parallel computing.

``` r
library(GA)
library(dplyr)
source("../scripts/GA_optimization_helper_functions.R")

#load data
COVID_contrasts = readRDS('../../data/optimization_datasets/COVID_contrasts.RDS')
COVID_contrasts_dict = readRDS('../../data/optimization_datasets/COVID_contrasts_dict.RDS')

non_COVID_contrasts = readRDS('../../data/optimization_datasets/all_curated_contrasts.RDS')
non_COVID_contrasts_dict = readRDS('../../data/optimization_datasets/non_COVID_contrasts_dict.RDS')

#select discovery data
COVID_contrasts = COVID_contrasts[which(COVID_contrasts_dict$use == 'discovery')]
COVID_contrasts_dict = COVID_contrasts_dict %>% dplyr::filter(use == 'discovery')

non_COVID_contrasts = non_COVID_contrasts[which(non_COVID_contrasts_dict$use == 'discovery')]
non_COVID_contrasts_dict = non_COVID_contrasts_dict %>% dplyr::filter(use == 'discovery')

#load prior info
prior_info_matrix = readRDS('../../data/optimization_datasets/prior_info_matrix.RDS')
filter_object = readRDS('../../data/optimization_datasets/filter_object.RDS')
pool_of_genes = c(filter_object$posGeneNames, filter_object$negGeneNames)

#initialize population of solutions
population_size = 10
expected_signature_size = 20
initial_population = initialize_population(population_size, solution_size = length(pool_of_genes), expected_signature_size, filter_object)

#generate all convex combinations of the weights
weight_grid = create_weight_grid(from = 0, to = 1, by = 0.5)

#run GA on the grid of weights 
maxiter = 5
GA_optimization_results = apply(weight_grid[1:2, ], 1,
                                function(x)
                                  ga(type = "binary",
                                     fitness = compute_fitness_from_binary,
                                     pool_of_genes = pool_of_genes,
                                     filter_object = filter_object,
                                     weight_COVID_vs_healthy = x[1],
                                     weight_COVID_vs_infection = x[2],
                                     weight_Non.Microbial = x[3],
                                     weight_Other = x[4],
                                     weight_Resp = x[5],
                                     prior_info_matrix = prior_info_matrix,
                                     weight_prior_info = x[6],
                                     max_signature_size = expected_signature_size,
                                     nBits = length(pool_of_genes),
                                     popSize = nrow(initial_population),
                                     maxiter = maxiter,
                                     suggestions = initial_population,
                                     names = pool_of_genes))

# binary solutions
lapply(GA_optimization_results, function(x) x@solution)
```

    ## $x_1_0_0_0_0_0
    ##      SIGMAR1 JSRP1 TRIOBP CD82 PIF1 TUBB4B ITM2C CPNE5 TMEM106C SRPRB FOXM1
    ## [1,]       0     0      0    0    0      0     0     0        0     0     0
    ##      POU2AF1 AGBL5 PSMD2 DYNLRB1 ASNA1 FARSB RANBP1 YIF1B ANXA7 TUBA4A DAD1
    ## [1,]       0     0     0       0     0     0      0     0     0      0    0
    ##      PREB ZSWIM1 MARCH2 FDFT1 NELFE MZB1 PPP1R7 ARHGDIA STMN1 EIF2B2 LYPLA2
    ## [1,]    0      0      0     0     0    0      0       0     0      0      0
    ##      YIF1A YARS IRAK1 B4GALT2 DERL3 TRAPPC1 STX5 MRPL37 PTTG1 CLPP BHLHA15
    ## [1,]     0    0     0       0     0       0    0      0     0    1       0
    ##      TNFRSF13B CD70 PHGDH PMVK MRPL12 ECH1 MRPL4 NDUFB9 IGLL5 TONSL NPLOC4
    ## [1,]         0    0     1    0      0    0     0      0     0     0      0
    ##      MRPS18B FEN1 GPRC5D TGFB1I1 POLR2L CACNA1F NRGN NUTF2 HDGF SAMD14 INF2
    ## [1,]       0    0      0       0      0       0    0     0    0      0    0
    ##      MYL6B NINJ1 PAFAH1B3 HIST1H2BJ HIST1H1C SWI5 PEX6 ANAPC11 BRMS1 ALYREF
    ## [1,]     0     0        0         1        0    0    0       0     0      0
    ##      FBXO18 CDC20 ETFB AURKAIP1 PSMD1 DDX41 STOML2 CCDC69 CSTB GPANK1 FAM207A
    ## [1,]      0     0    0        0     0     0      0      0    0      0       0
    ##      AQP3 NT5M PRDX2 PTCRA HMGN2 CMTM5 FERMT3 NT5DC2 PDAP1 SUGP1 ALOX5AP TCEAL3
    ## [1,]    0    0     0     0     0     0      0      0     0     0       0      0
    ##      H1FX S100A10 DNAJC4 UQCRC1 BIRC5 CHMP6 PLK3 MYBL2 CHST13 HCFC1R1 CLEC2L
    ## [1,]    0       0      0      0     0     0    0     0      0       0      0
    ##      POP7 ALOX12 VSIG2 DLGAP5 ATP5G1 UCHL1 BANF1 TXNDC5 CCNB1 VMAC LMAN2 UBE2L3
    ## [1,]    0      0     0      0      0     0     0      0     0    0     0      0
    ##      SMPD4 UROD TMED1 MAOB ZBTB32 SLC39A3 RAD23A NDUFA3 NDUFA8 ADI1 XBP1 SZRD1
    ## [1,]     0    0     0    0      0       0      0      0      1    0    0     0
    ##      NDUFB10 SNRPC DDX49 RBM42 POLDIP2 SNF8 PGAM5 B4GALT3 KHSRP PF4 SPCS1
    ## [1,]       0     0     0     0       0    0     0       0     0   1     0
    ##      C1orf122 AQP10 SMIM1 CCDC124 COMMD1 TMED9 AIFM2 C21orf58 SPCS2 SMOX GTF2A2
    ## [1,]        0     0     0       0      0     0     0        0     0    1      0
    ##      NPDC1 CLTA C15orf52 NUDT18 TCEAL4 RUFY1 AP2M1 DPCD PTMS CENPO PSMD11 DCTN3
    ## [1,]     0    0        0      0      0     0     0    1    0     0      0     0
    ##      TACO1 NDUFB11 VKORC1 EHD3 HINT2 TREML1 CUEDC2 CRAT MANF DHFR ORMDL2 EMP3
    ## [1,]     0       0      0    1     0      0      0    0    0    0      0    0
    ##      PLOD1 HMBS LRRC32 DOHH SH3GL1 PES1 ALG1 PSMD3 DNM2 SIRT3 PAICS LSM2
    ## [1,]     0    0      0    0      0    0    0     0    0     0     0    0
    ##      HNRNPAB RITA1 EMC4 RFC2 YIPF3 SLC7A5 MRPS12 MCM6 CHID1 CHAF1A EIF4A3 SGTA
    ## [1,]       0     0    0    0     0      0      0    0     0      0      0    0
    ##      YIPF2 MYO1D RRP1 GUCD1 FAXDC2 MYL9 CHCHD3 ATAD3A PPIB TUBA1C ITPR2 TAOK1
    ## [1,]     0     0    0     1      0    0      0      0    0      0     0     0
    ##      MYCBP2 LNPEP MS4A7 BIRC6 CREB1 PCGF3 GPATCH2L PKN2 ACAP2 PUM1 BRWD3 TAF1
    ## [1,]      0     0     0     0     0     0        0    0     0    0     0    0
    ##      RC3H1 COL4A3BP GSK3B UBR2 ZFX ATL3 DDX3X DPP8 BRAF PRRC2C GAPVD1 CREBRF
    ## [1,]     0        0     0    1   0    0     0    0    0      0      0      1
    ##      KPNA4 PIKFYVE PHC3 ZNF217 BOD1L1 FAM126B USP9X TTN DOCK5 KMT2C ASH1L
    ## [1,]     0       0    0      0      0       0     0   0     0     0     0
    ##      TRIM33 BTAF1 ZC3H7A CEP192 PRKAR1A DCP2 ANKRD44 UBR5 GPR155 ZFC3H1 DCAF4L1
    ## [1,]      0     0      0      0       0    0       0    0      0      0       0
    ##      IPCEF1 LPGAT1 AGO3 RICTOR TRIM44 TMEM154 CD44 MACF1 ATAD2B SYNE2 GSAP
    ## [1,]      0      0    0      0      0       0    0     0      0     0    1
    ##      MIER1 MAP3K2 UBN2 SP4 FAM26F SCAF11 GBP2 AHSA2 MARCH6 ASXL2 TMEM65 ARFGEF1
    ## [1,]     0      0    0   1      0      0    0     0      0     0      0       0
    ##      C16orf72 ARID2 HERC4 TNRC6B CREBZF BAZ2B RSBN1L AQR ZBTB18 MSL1 SMAD4
    ## [1,]        0     0     0      0      0     0      0   0      0    0     0
    ##      TMOD2 NIPBL MLLT10 SMAD2 SIKE1 ZNF431 HMBOX1 GOLGA6L9 CCNG2 FAM111A HIVEP2
    ## [1,]     0     0      0     0     0      0      0        0     0       0      1
    ##      SMG1 KMT2A AKAP9 MSL2 ZNF117 SMIM14 NR3C1 SETD5 PAN3 CNOT6L CCSAP P2RY14
    ## [1,]    0     0     0    1      0      0     0     0    0      0     0      0
    ##      ZMYM6 NXPE3 SCLT1 MED23 WDR11 ZNF493 CARD8 SPAG9 ATP2B4 SORL1 RIF1 CEP135
    ## [1,]     0     0     0     0     0      0     0     0      1     0    0      0
    ##      TGFBR1 NKTR GTF2A1 SUDS3 GSE1 ANGPT2 SON VEZF1 CTDSPL2 RGS17 VPS41 VCPIP1
    ## [1,]      0    0      0     0    0      0   0     0       0     0     0      0
    ##      ELF2 TVP23B YTHDC2 MED28 GPR65 CARF SLK FOXK1 ANKRD17 RAB11FIP2 SLC25A46
    ## [1,]    0      0      0     0     0    0   0     0       0         0        0
    ##      CLOCK ZNF264 MRAS KATNBL1 UBE2W ARID1A HEG1 DCP1A ZCCHC11 ZNF91 PDP1
    ## [1,]     0      0    1       0     0      0    0     0       0     0    0
    ##      ZNF770 DR1 TMEM170A TTC39B KPNB1 ARAP2 NR2C2 GAN MDN1 SLC26A2 SLC35A3
    ## [1,]      0   0        0      1     0     0     0   0    0       0       0
    ##      USP53 PHIP ZDHHC17 ATM DIP2A JMY CFLAR MBNL1 PHF3 ROCK2 RC3H2 TRPM7 ATF7IP
    ## [1,]     0    0       0   0     0   1     0     0    0     0     0     0      0
    ##      MED13 PLEKHA2 TTC14 ZC3H11A N4BP2L2 PCF11 FRYL FAM8A1 TRIM56
    ## [1,]     0       0     0       0       0     0    0      0      0
    ## 
    ## $x_0_1_0_0_0_0
    ##      SIGMAR1 JSRP1 TRIOBP CD82 PIF1 TUBB4B ITM2C CPNE5 TMEM106C SRPRB FOXM1
    ## [1,]       0     0      0    0    0      0     0     0        0     0     0
    ##      POU2AF1 AGBL5 PSMD2 DYNLRB1 ASNA1 FARSB RANBP1 YIF1B ANXA7 TUBA4A DAD1
    ## [1,]       0     0     0       0     0     0      0     0     0      0    0
    ##      PREB ZSWIM1 MARCH2 FDFT1 NELFE MZB1 PPP1R7 ARHGDIA STMN1 EIF2B2 LYPLA2
    ## [1,]    0      0      0     0     0    0      0       0     0      0      0
    ##      YIF1A YARS IRAK1 B4GALT2 DERL3 TRAPPC1 STX5 MRPL37 PTTG1 CLPP BHLHA15
    ## [1,]     0    0     0       0     0       0    0      0     0    0       0
    ##      TNFRSF13B CD70 PHGDH PMVK MRPL12 ECH1 MRPL4 NDUFB9 IGLL5 TONSL NPLOC4
    ## [1,]         0    0     0    0      0    0     0      0     0     0      1
    ##      MRPS18B FEN1 GPRC5D TGFB1I1 POLR2L CACNA1F NRGN NUTF2 HDGF SAMD14 INF2
    ## [1,]       0    0      0       0      0       0    0     0    1      0    0
    ##      MYL6B NINJ1 PAFAH1B3 HIST1H2BJ HIST1H1C SWI5 PEX6 ANAPC11 BRMS1 ALYREF
    ## [1,]     0     0        0         0        0    0    0       0     0      0
    ##      FBXO18 CDC20 ETFB AURKAIP1 PSMD1 DDX41 STOML2 CCDC69 CSTB GPANK1 FAM207A
    ## [1,]      0     0    0        0     0     0      0      0    0      0       0
    ##      AQP3 NT5M PRDX2 PTCRA HMGN2 CMTM5 FERMT3 NT5DC2 PDAP1 SUGP1 ALOX5AP TCEAL3
    ## [1,]    0    0     0     0     0     0      1      0     0     0       1      0
    ##      H1FX S100A10 DNAJC4 UQCRC1 BIRC5 CHMP6 PLK3 MYBL2 CHST13 HCFC1R1 CLEC2L
    ## [1,]    0       0      0      0     0     0    0     0      0       0      0
    ##      POP7 ALOX12 VSIG2 DLGAP5 ATP5G1 UCHL1 BANF1 TXNDC5 CCNB1 VMAC LMAN2 UBE2L3
    ## [1,]    0      0     0      0      0     0     0      0     0    0     0      0
    ##      SMPD4 UROD TMED1 MAOB ZBTB32 SLC39A3 RAD23A NDUFA3 NDUFA8 ADI1 XBP1 SZRD1
    ## [1,]     0    0     0    0      1       0      0      0      0    0    0     0
    ##      NDUFB10 SNRPC DDX49 RBM42 POLDIP2 SNF8 PGAM5 B4GALT3 KHSRP PF4 SPCS1
    ## [1,]       0     0     0     0       0    0     0       0     0   0     0
    ##      C1orf122 AQP10 SMIM1 CCDC124 COMMD1 TMED9 AIFM2 C21orf58 SPCS2 SMOX GTF2A2
    ## [1,]        0     0     0       0      0     1     0        0     1    0      0
    ##      NPDC1 CLTA C15orf52 NUDT18 TCEAL4 RUFY1 AP2M1 DPCD PTMS CENPO PSMD11 DCTN3
    ## [1,]     0    0        0      0      0     0     0    0    0     0      0     0
    ##      TACO1 NDUFB11 VKORC1 EHD3 HINT2 TREML1 CUEDC2 CRAT MANF DHFR ORMDL2 EMP3
    ## [1,]     0       0      0    0     0      0      0    0    0    0      0    0
    ##      PLOD1 HMBS LRRC32 DOHH SH3GL1 PES1 ALG1 PSMD3 DNM2 SIRT3 PAICS LSM2
    ## [1,]     0    0      0    0      0    0    0     0    0     0     0    0
    ##      HNRNPAB RITA1 EMC4 RFC2 YIPF3 SLC7A5 MRPS12 MCM6 CHID1 CHAF1A EIF4A3 SGTA
    ## [1,]       0     1    0    0     0      0      0    0     0      0      0    0
    ##      YIPF2 MYO1D RRP1 GUCD1 FAXDC2 MYL9 CHCHD3 ATAD3A PPIB TUBA1C ITPR2 TAOK1
    ## [1,]     0     0    0     0      1    0      0      1    0      0     0     0
    ##      MYCBP2 LNPEP MS4A7 BIRC6 CREB1 PCGF3 GPATCH2L PKN2 ACAP2 PUM1 BRWD3 TAF1
    ## [1,]      0     0     0     0     0     1        0    0     0    0     0    0
    ##      RC3H1 COL4A3BP GSK3B UBR2 ZFX ATL3 DDX3X DPP8 BRAF PRRC2C GAPVD1 CREBRF
    ## [1,]     0        0     0    0   0    0     0    0    0      0      0      0
    ##      KPNA4 PIKFYVE PHC3 ZNF217 BOD1L1 FAM126B USP9X TTN DOCK5 KMT2C ASH1L
    ## [1,]     0       0    0      0      0       0     0   0     0     0     0
    ##      TRIM33 BTAF1 ZC3H7A CEP192 PRKAR1A DCP2 ANKRD44 UBR5 GPR155 ZFC3H1 DCAF4L1
    ## [1,]      0     0      0      0       0    0       0    0      0      0       0
    ##      IPCEF1 LPGAT1 AGO3 RICTOR TRIM44 TMEM154 CD44 MACF1 ATAD2B SYNE2 GSAP
    ## [1,]      0      0    0      0      0       1    0     0      0     0    0
    ##      MIER1 MAP3K2 UBN2 SP4 FAM26F SCAF11 GBP2 AHSA2 MARCH6 ASXL2 TMEM65 ARFGEF1
    ## [1,]     0      0    0   0      0      0    0     0      0     0      0       0
    ##      C16orf72 ARID2 HERC4 TNRC6B CREBZF BAZ2B RSBN1L AQR ZBTB18 MSL1 SMAD4
    ## [1,]        0     0     0      0      0     1      0   0      0    0     0
    ##      TMOD2 NIPBL MLLT10 SMAD2 SIKE1 ZNF431 HMBOX1 GOLGA6L9 CCNG2 FAM111A HIVEP2
    ## [1,]     0     0      0     0     0      0      0        0     0       0      0
    ##      SMG1 KMT2A AKAP9 MSL2 ZNF117 SMIM14 NR3C1 SETD5 PAN3 CNOT6L CCSAP P2RY14
    ## [1,]    1     0     0    0      0      0     0     0    0      0     0      1
    ##      ZMYM6 NXPE3 SCLT1 MED23 WDR11 ZNF493 CARD8 SPAG9 ATP2B4 SORL1 RIF1 CEP135
    ## [1,]     0     0     0     0     1      0     0     0      0     0    0      0
    ##      TGFBR1 NKTR GTF2A1 SUDS3 GSE1 ANGPT2 SON VEZF1 CTDSPL2 RGS17 VPS41 VCPIP1
    ## [1,]      0    0      0     0    0      0   0     0       1     0     0      0
    ##      ELF2 TVP23B YTHDC2 MED28 GPR65 CARF SLK FOXK1 ANKRD17 RAB11FIP2 SLC25A46
    ## [1,]    1      0      0     0     0    0   0     0       0         1        0
    ##      CLOCK ZNF264 MRAS KATNBL1 UBE2W ARID1A HEG1 DCP1A ZCCHC11 ZNF91 PDP1
    ## [1,]     0      0    0       0     0      0    0     0       0     0    0
    ##      ZNF770 DR1 TMEM170A TTC39B KPNB1 ARAP2 NR2C2 GAN MDN1 SLC26A2 SLC35A3
    ## [1,]      0   0        0      0     0     0     0   0    0       0       0
    ##      USP53 PHIP ZDHHC17 ATM DIP2A JMY CFLAR MBNL1 PHF3 ROCK2 RC3H2 TRPM7 ATF7IP
    ## [1,]     0    0       0   0     0   0     0     0    0     1     0     0      0
    ##      MED13 PLEKHA2 TTC14 ZC3H11A N4BP2L2 PCF11 FRYL FAM8A1 TRIM56
    ## [1,]     0       0     0       0       0     0    0      0      0
