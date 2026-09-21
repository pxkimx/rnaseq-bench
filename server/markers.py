"""Reference marker panels for cluster annotation and cell-cycle scoring."""

CELL_TYPE_MARKERS = {
    # immune (PBMC / blood / tumour infiltrate)
    "CD4+ T": ["IL7R", "CCR7", "CD3E", "CD3D", "CD4", "TCF7", "LEF1", "MAL"],
    "CD8+ T": ["CD8A", "CD8B", "CD3E", "GZMK", "CCL5", "GZMH"],
    "Regulatory T": ["FOXP3", "IL2RA", "CTLA4", "TNFRSF4"],
    "γδ / MAIT T": ["TRDC", "TRGC1", "KLRB1", "SLC4A10"],
    "NK": ["GNLY", "NKG7", "PRF1", "KLRD1", "GZMB", "KLRF1", "NCAM1"],
    "B": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "PAX5"],
    "Plasma": ["JCHAIN", "MZB1", "XBP1", "IGHG1", "SDC1"],
    "CD14+ Monocyte": ["CD14", "LYZ", "S100A8", "S100A9", "FCN1", "VCAN"],
    "FCGR3A+ Monocyte": ["FCGR3A", "MS4A7", "LST1", "CDKN1C", "IFITM3"],
    "cDC": ["FCER1A", "CST3", "CLEC10A", "CD1C", "CLEC9A"],
    "pDC": ["LILRA4", "IL3RA", "CLEC4C", "TCF4"],
    "Macrophage": ["C1QA", "C1QB", "C1QC", "CD68", "APOE", "MRC1"],
    "Neutrophil": ["CSF3R", "FCGR3B", "CXCL8", "S100A12"],
    "Mast": ["TPSAB1", "CPA3", "KIT", "MS4A2"],
    "Platelet": ["PPBP", "PF4", "GP9", "TUBB1"],
    "Erythroid": ["HBB", "HBA1", "HBA2", "ALAS2"],
    "HSPC": ["CD34", "PROM1", "SPINK2"],
    # tissue
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "KDR"],
    "Smooth muscle / pericyte": ["ACTA2", "MYH11", "TAGLN", "RGS5"],
    "Hepatocyte": ["ALB", "APOA1", "TTR", "SERPINA1"],
    "Keratinocyte": ["KRT5", "KRT14", "KRT1", "KRT10"],
    "Melanocyte": ["PMEL", "MLANA", "TYR"],
    "Neuron": ["SNAP25", "SYT1", "RBFOX3", "STMN2"],
    # brain barriers: the choroid plexus and ependyma are epithelia, but they share almost nothing with
    # the EPCAM/keratin panel above, so without these they come back as "Unassigned"
    "Choroid plexus epithelium": ["TTR", "FOLR1", "CLIC6", "HTR2C", "KCNJ13", "OTX2", "PRLR", "AQP1"],
    "Ependymal": ["FOXJ1", "PIFO", "CCDC153", "TMEM212", "HDC", "SPAG17"],
    # pericytes are usually reported separately from smooth muscle; these markers exclude the
    # contractile genes (ACTA2/MYH11/TAGLN) that define a true VSMC
    "Pericyte": ["PDGFRB", "KCNJ8", "HIGD1B", "ABCC9", "PTH1R", "RGS5"],
    "Astrocyte": ["GFAP", "AQP4", "SLC1A3", "ALDH1L1"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MOBP"],
    "OPC": ["PDGFRA", "OLIG1", "OLIG2", "CSPG4"],
    "Microglia": ["P2RY12", "CX3CR1", "TMEM119", "C1QA"],
    "Cardiomyocyte": ["TNNT2", "MYH6", "MYL7", "ACTC1"],
    # STMN1 and TUBA1B used to be in here and should not have been: tubulin and stathmin are broadly
    # expressed (TUBA1B in 40-90% of cells in a post-mitotic epithelium), so they dominated the score and
    # labelled ordinary clusters as cycling. These are specific to cells actually in S/G2/M.
    "Proliferating": ["MKI67", "TOP2A", "CENPF", "UBE2C", "BIRC5", "CCNB1", "CDK1", "ASPM"],
    # vascular / aorta
    "Modulated SMC (fibromyocyte)": ["LGALS3", "VCAM1", "FN1", "TNFRSF11B", "SPP1", "LUM"],
    "Adventitial fibroblast": ["PI16", "DPT", "CD34", "PDGFRA", "GSN"],
    "Mesothelial": ["MSLN", "UPK3B", "WT1", "KRT19"],
    "Adipocyte": ["ADIPOQ", "LEP", "PLIN1", "CFD"],
    "Schwann / glia": ["PLP1", "SOX10", "MPZ", "S100B"],
    "Lymphatic endothelial": ["PROX1", "LYVE1", "FLT4", "CCL21"],
}

# Tirosh et al. 2016 (as used in the Seurat / Scanpy cell-cycle tutorials)
S_GENES = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7",
           "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN",
           "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45",
           "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B",
           "BRIP1", "E2F8"]
G2M_GENES = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2",
             "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2",
             "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1",
             "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
             "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
             "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"]


# Dissociation-induced stress response (van den Brink et al., Nat Methods 2017). Warm protease
# dissociation makes cells transcribe immediate-early and heat-shock genes within minutes; the result
# looks like a distinct population and is routinely mistaken for one. Jovic et al. (Clin Transl Med
# 2022) name this as a principal pitfall of scRNA-seq and recommend dissociating at 4 degrees.
# Deliberately limited to immediate-early and chaperone genes: the fuller published list includes
# genes such as DCN, DES and MT1/MT2 that are genuine markers of fibroblasts, muscle and other
# populations, and scoring those would flag real cell types as stressed.
DISSOCIATION_GENES = [
    "FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "ATF3", "IER2", "IER3", "IER5",
    "DUSP1", "DUSP5", "ZFP36", "ZFP36L1", "SOCS3", "KLF2", "KLF4", "KLF6", "BTG1", "BTG2", "NR4A1",
    "PPP1R15A", "CEBPB", "CEBPD", "RHOB", "SGK1", "MCL1", "NFKBIA", "CCN1", "CYR61", "SQSTM1",
    "HSPA1A", "HSPA1B", "HSPA8", "HSPB1", "HSPH1", "HSP90AA1", "HSP90AB1", "DNAJA1", "DNAJB1",
    "BAG3", "UBC", "MYC", "PER1", "TNFAIP3",
]
