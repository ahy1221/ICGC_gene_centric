PCAWG.driver <- new.env(parent = .GlobalEnv)
TCGA.driver <- new.env(parent = .GlobalEnv)
add_geneset <- function(genes, env, name) {
  env.name = deparse(substitute(env))
  INFO = sprintf("Add %s to %s, including %i genes...", name, env.name, length(genes))
  print(INFO)
  assign(name, genes, envir = env)
}
add_geneset(TCGA.HCC, TCGA.driver, "HCC")
liver.driver.TCGA <- c("TP53", "CTNNB1", "ALB", "AXIN1", "BAP1", "KEAP1", "NFE2L2", "LZTR1", "RB1",
                       "PIK3CA", "AZIN1", "KRAS", "IL6ST", "CDKN2A", "EEF1A1",
                       "ARID2", "ARID1A", "GPATCH4", "ACVR2A", "APOB", "CREB3L3", "NRAS", "AHCTF1",
                       "HIST1H1C")




PCAWGA.drivers <-c("
ALB
APC
ARID1A
ARID2
AXIN1
CDKN2A
CTNNB1
KRAS
PBRM1
PIK3CA
PTEN
RB1
SETD2
SMAD4
TP53
VHL
RNF43
MAP2K4
PIK3R1
CBFB
BAP1
ARHGAP35
SPOP
BRAF
CDH1
FBXW7
FOXA1
RBM10
RPL22
ACVR2A
NRAS
MAP3K1
KMT2C
SF3B1
CDK12
ATM
ZFP36L2
NFE2L2
KEAP1
CDKN1A
AKT1
BRD7
GATA3
ACVR1B
SOX9
B2M
MAP2K7
KDM6A
TGFBR2
EPHA2
TBX3
CDKN1B
CASP8
TBL1XR1
SETDB1
RPS6KA3
TP53
SMAD4
KRAS
ARID1A
PTEN
TP53
KDM6A
CDKN1A
FGFR3
ARID1A
RB1
TP53
TP53
RB1
CBFB
GATA3
MAP3K1
PIK3CA
TP53
PTEN
AKT1
CDH1
MAP2K4
CTCF
SMAD4
CBFB
CDH1
GATA3
MAP3K1
PIK3CA
TP53
PTEN
AKT1
KMT2C
MAP2K4
PIK3R1
FOXA1
NOTCH2
NCOR1
CTCF
SMAD4
EGFR
TP53
PTEN
IDH1
RB1
DDX3X
PTCH1
CTDNEP1
KMT2D
SMO
PRKAR1A
LDB1
TCF4
KBTBD4
KDM6A
KMT2C
PIK3CA
SMARCA4
BRPF1
FBXW7
TBR1
BRAF
FGFR1
NF1
DDX3X
EGFR
IDH1
PTCH1
TP53
PTEN
CIC
CTDNEP1
BRAF
PIK3CA
KMT2D
SMO
ATRX
LDB1
FGFR1
KBTBD4
BCOR
PRKAR1A
FUBP1
PTPN11
NF1
KMT2C
KDM6A
ALB
APC
ARID1A
AXIN1
CDKN2A
CTNNB1
FBXW7
KRAS
MEN1
NFE2L2
PBRM1
PIK3CA
PTEN
RB1
SETD2
SMAD4
TP53
VHL
ARID2
CDH1
BAP1
CBFB
PIK3R1
KMT2D
MAP2K4
KDM6A
RNF43
CASP8
CDKN1A
SPOP
BRAF
RBM10
KMT2C
ATM
GATA3
ARHGAP35
MAP3K1
RPL22
NRAS
ACVR2A
HLA-A
KEAP1
B2M
ZFP36L2
ACVR1B
SF3B1
NOTCH1
STK11
HRAS
FOXA1
AKT1
BRD7
ELF3
DAXX
TGFBR2
CDKN1B
EPHA2
SOX9
BAX
NF1
MAP2K7
TSC1
RPS6KA3
CTCF
SMARCB1
APC
KRAS
TP53
PIK3CA
SOX9
SMAD4
FBXW7
PCBP1
CTNNB1
ZFP36L2
SMAD2
B2M
TCF7L2
ACVR2A
NRAS
PTEN
ALB
APC
ARID1A
ARID2
AXIN1
CDKN2A
CTNNB1
KRAS
PIK3CA
PTEN
SMAD4
TP53
RNF43
RPL22
MAP2K4
ZFP36L2
ACVR2A
RB1
SOX9
NFE2L2
BAP1
RPS6KA3
CDKN1A
B2M
BRD7
RBM10
NRAS
BRAF
PBRM1
FBXW7
TGFBR2
ATM
HLA-A
KDM6A
MAP2K7
HIST1H4D
ACVR1B
ARID1A
CDKN2A
SMAD4
TP53
KRAS
APC
PIK3CA
PTEN
MAP2K7
CBFB
GATA3
MAP3K1
PIK3CA
PIK3R1
PTEN
TP53
CDH1
AKT1
RB1
ARID1A
KRAS
FBXW7
PPP2R1A
ARHGAP35
MAP2K4
CTCF
KMT2C
NF1
ATRX
EGFR
IDH1
TP53
PTEN
BRAF
CIC
FGFR1
NF1
PTPN11
FUBP1
PIK3CA
TP53
CDKN2A
PIK3CA
AJUBA
NOTCH1
CASP8
HRAS
PA2G4
B2M
BTG1
CCND3
CREBBP
EZH2
HIST1H1C
KMT2D
MYD88
TNFRSF14
TP53
TMSB4X
HIST1H2BC
GNAI2
CARD11
KLHL6
STAT6
DDX3X
ARID1A
GRB2
FBXO11
MEF2B
NFKBIE
RRAGC
GNA13
MCL1
PAX5
SMARCA4
STAT3
TMEM30A
PRKCD
RFTN1
HIST1H2AG
HIST1H1B
RHOA
TBL1XR1
HLA-B
NRAS
SF3B1
SRSF7
PTEN
PHF6
ITPKB
NXF1
TET2
ATM
EBF1
TP53
PBRM1
SETD2
VHL
BAP1
PTEN
KDM5C
MTOR
PBRM1
SETD2
VHL
BAP1
PTEN
TP53
MTOR
ALB
ARID2
AXIN1
CTNNB1
TP53
ARID1A
NFE2L2
RB1
CDKN1A
BRD7
BAP1
RPS6KA3
KEAP1
PTEN
RPL22
APOB
DYRK1A
HNF4A
ACVR2A
EEF1A1
SETDB1
CAMK1
RPL5
TP53
KRAS
STK11
KEAP1
RBM10
NFE2L2
TP53
CDKN2A
NOTCH1
KMT2D
CREBBP
PIK3CA
NFE2L2
TP53
CDKN2A
KRAS
STK11
NOTCH1
KEAP1
RBM10
PIK3CA
CTC-512J12.6
B2M
BTG1
CCND3
CREBBP
HIST1H1C
KMT2D
TNFRSF14
TP53
GNAI2
HIST1H2BC
EZH2
TMSB4X
MYD88
STAT6
FBXO11
MEF2B
CARD11
RRAGC
ARID1A
GRB2
KLHL6
DDX3X
GNA13
SMARCA4
HIST1H1B
HLA-B
PRKCD
RHOA
RFTN1
STAT3
TMEM30A
TBL1XR1
SRSF7
PTEN
NFKBIE
EBF1
SF3B1
TP53
NXF1
ATM
MYD88
NFKBIE
BRAF
ZNF292
IRF4
B2M
BTG1
CCND3
CREBBP
HIST1H1C
KMT2D
MYD88
TNFRSF14
TP53
GNAI2
HIST1H2BC
TMSB4X
KLHL6
EZH2
CARD11
STAT6
DDX3X
NFKBIE
ARID1A
MEF2B
GRB2
FBXO11
RRAGC
GNA13
MCL1
PAX5
TMEM30A
SMARCA4
STAT3
RFTN1
TBL1XR1
HIST1H1B
RHOA
HLA-B
PRKCD
SF3B1
SRSF7
PTEN
HIST1H2AG
ITPKB
NXF1
ATM
EBF1
TET2
DNMT3A
EZH2
NRAS
IDH2
DNMT3A
TP53
RB1
CDK12
LATS1
APC
ARID1A
CDKN2A
CTNNB1
IDH1
KRAS
NFE2L2
PBRM1
PIK3CA
PTEN
RB1
SMAD4
TP53
VHL
SETD2
CDH1
BRAF
FBXW7
MEN1
BAP1
MAP2K4
ARID2
KDM6A
B2M
RPL22
CDKN1A
KEAP1
ZFP36L2
ALB
RNF43
ELF3
NOTCH1
ACVR1B
NRAS
CBFB
SPOP
SMARCA4
MAP3K1
ACVR2A
PIK3R1
ARHGAP35
ATM
AXIN1
CASP8
DDX3X
HLA-A
KMT2C
KMT2D
RBM10
TGFBR2
SF3B1
EPHA2
PTCH1
MAP2K7
GATA3
CDKN1B
NF1
FOXA1
U2AF1
HRAS
RPL5
KDM5C
EGFR
ARID1A
CDKN2A
KRAS
SMAD4
TP53
ZFP36L2
MAP2K4
RNF43
TGFBR2
GNAS
ACVR1B
RBM10
RIPK4
KDM6A
SETD2
ACVR2A
DAXX
MEN1
PTEN
ATRX
SETD2
DYNC1I1
SPOP
TP53
PTEN
FOXA1
TP53
RB1
H3F3B
H3F3A
BRAF
CDKN2A
NRAS
TP53
PTEN
RAC1
NF1
CDKN2A
NFE2L2
NOTCH1
PIK3CA
TP53
KMT2D
HLA-A
AJUBA
FAT1
CASP8
HRAS
KRAS
CTC-512J12.6
TP53
PIK3CA
ARID1A
SMAD4
CTNNB1
PLK1
RELA
PTEN
BRAF
PIK3CA
PIK3R1
PTEN
TP53
PPP2R1A
ARHGAP35
ARID1A
CTNNB1
KRAS
FBXW7")
PCAWG.drivers = unlist(strsplit(PCAWGA.drivers, "\n"))

          