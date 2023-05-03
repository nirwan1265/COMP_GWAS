# Interpretation sets: #######

# Soluble phosphorus
# add units and citation

feat_sol <- c(
  "lab",
  # Labile Pi (Yang2014)
  
  "PBR1",
  # The amount of phosphorous using the Bray1 method  (Shangguan2014)
  # (more successful in acid soils)
  
  "PHO1", 
  # The amount of water soluble phosphorous (Shangguan2014)
  
  "ext"
  #  Extractable phosphorus (Miller2021)
)

# cols.cor[feat_sol,feat_sol]
# cor_avg(cols.cor,feat_sol)

# Phosphorus retention,
# check scaling for correct units

feat_ret <- c(
  "occ",  
  # Occluded P (Yang2014)
  # occluded phosphorus usually in clays - high Al soils 
  
  "ret_VH",    
  # Probability of Very High Retention soil (Batjes2011)
  # probability of finding a soil type/class with expected Very High P retention
  # This class has higher proportions of Ferralsols - Andosols 
  
  "PNZ1"  
  # Phosphorous retention by New Zealand method (Shangguan2014)
  # as weight percentage
) 

# cols.cor[feat_ret,feat_ret]
# cor_avg(cols.cor,feat_ret)

# Total phosphorus 
# ppm weight

feat_tot <- c(
  # TotalPi (Yang2014)
  "tot",
  # Total phosphorus (Shangguan2014)      
  "TP1",
  # Total phosphorus (He2021)
  "stp10"
)
# cols.cor[feat_tot,feat_tot]
# cor_avg(cols.cor,feat_tot)


# Alkaline soil phosphorus
feat_alk <- c(
  "POL1",   
  # The amount of phosphorous by Olsen method (Shangguan2014)
  # Olsen method works better on alkali soils
  
  "ret_Mo"  
  # This soil class has higher proportion of alkali soils
) 

# cols.cor[feat_alk,feat_alk]
# cor_avg(cols.cor,feat_alk)

# P/N ratios

feat_PNrat <- c(
  "lab_N",
  
  "PBR1_N",
  
  "PHO1_N",
  
  "ext_N",
  
  "NPlim"
  # Nitrogen vs phosphorus limitation (Du2020)
  # ln-transformed NREdom/PREdom
  # not a ratio but zn estimate of whether Nitrogen is more likely to be limiting
  # < -0.16 to indicate P limitation
  # > 0.16 to indicate N limitation
)  

