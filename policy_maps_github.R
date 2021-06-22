#################
# Maps for health policy and IMR
# Citation: Webster JL, Paul D, Purtle J, Locke R, Goldstein ND. State-level Social and Economic Policies and their Association with Perinatal and Infant Outcomes. Manuscript in preparation.
# 11/25/20 -- Neal Goldstein
#################


### FUNCTIONS ###

library(psych) #describe, describeBy
library(gmodels) #CrossTable
library("tidyLPA") #latent profile analysis
library("rgdal") #readOGR
library("sp") #shapefile
library("maptools") #manipulation functions
library("plotrix") #plotting functions
library("RColorBrewer") #color palette


### READ DATA ###

#Infant mortality from CDC Wonder (via NCHS, 2018): https://wonder.cdc.gov
##Select Infant Deaths (Linked Birth/Infant Death Records) > Linked Birth / Infant Death Records for 2017-2018 with ICD 10 codes (expanded)
##Group results by Maternal State, Additional Rates per 1,000, Infant characertistics>Year of Death to 2018, Precision to 1 decimal
#Fetal mortality from CDC Wonder (via NCHS, 2018): https://wonder.cdc.gov
##Select Fetal Deaths > Fetal Deaths for 2014 - 2018 (expanded)
##Group results by Maternal State, Delivery characteristics>Year to 2018, Precision to 1 decimal
#Live births from CDC Wonder (via NCHS, 2018): https://wonder.cdc.gov
##Select Births > Natality for 2016 - 2019 (expanded)
##Group results by Maternal State, Delivery characteristics>Year to 2018, Precision to 0 decimal
#Cigarette tax from Campaign for Tobacco-Free Kids (2020): https://www.tobaccofreekids.org/assets/factsheets/0097.pdf
#Minimum wage from US Dept of Labor (2020): https://www.dol.gov/agencies/whd/mw-consolidated and https://www.dol.gov/agencies/whd/minimum-wage/state
#Paid parental leave, score and grade, from the National Partnership for Women & Families (2005): https://www.leg.state.nv.us/App/NELIS/REL/79th2017/ExhibitDocument/OpenExhibitDocument?exhibitId=29512&fileDownloadName=0330ab266_ParentalLeaveReportMay05.pdf
#Paid parental leave, type, from the National Conference of State Legislators (2020): https://www.ncsl.org/Portals/1/Documents/cyf/Paid-Family-Leave_v05.pdf
#Tax credit from Tax Credits for Workers and Their Families (2019): http://www.taxcreditsforworkersandfamilies.org/state-tax-credits/
policy_data = read.csv("MCH_policy_data.csv", as.is=T, stringsAsFactors=F)

#state cartographic boundaries: https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html
us_carto = readOGR("cb_2018_us_state_20m/", "cb_2018_us_state_20m")


### EXPLORATORY FACTOR ANALYSIS ###

# #FA recodings
# policy_data$Wage_fed_comparison_7.25 = as.integer(unclass(as.factor(policy_data$Wage_fed_comparison_7.25)))
# #policy_data$Parental_leave_grade = as.integer(unclass(as.factor(policy_data$Parental_leave_grade)))
# policy_data$Parental_leave_type = as.integer(unclass(as.factor(policy_data$Parental_leave_type)))
# policy_data$Tax_EITC = as.integer(unclass(as.factor(policy_data$Tax_EITC)))
# policy_data$Tax_CTC = as.integer(unclass(as.factor(policy_data$Tax_CTC)))
# policy_data$Tax_CDCTC = as.integer(unclass(as.factor(policy_data$Tax_CDCTC)))
# 
# #FA variable list
# factor_vars = c("Cigarette_tax","Cigarette_rank","Wage","Wage_fed_comparison_7.25","Parental_leave_type","Tax_EITC","Tax_EITC_fedpercent","Tax_CTC","Tax_CDCTC","Tax_CDCTC_fedpercent")
# 
# factors = factanal(policy_data[,factor_vars], 3, scores=c("regression"), rotation="varimax")
# factors


### AND LATENT PROFILE ANALYSIS ###

#recodings
policy_data$Parental_leave_type_cat = ifelse(policy_data$Parental_leave_type=="Broad", 3, ifelse(policy_data$Parental_leave_type=="Narrow", 2, ifelse(policy_data$Parental_leave_type=="Pending", 1, 0)))
policy_data$Tax_cumulative_cat = (policy_data$Tax_EITC=="Yes") + (policy_data$Tax_CTC=="Yes") + (policy_data$Tax_CDCTC=="Yes")

#add factors to data
#policy_data = cbind(policy_data, factors$scores)

#scaling parameters for LPA plot
#policy_data$Cigarette_tax = policy_data$Cigarette_tax * 3
#policy_data$Parental_leave_type_cat = policy_data$Parental_leave_type_cat * 5
#policy_data$Tax_cumulative_cat = policy_data$Tax_cumulative_cat * 5

#LPA: https://cran.r-project.org/web/packages/tidyLPA/vignettes/Introduction_to_tidyLPA.html
#lpa = estimate_profiles(policy_data[, c("Factor1", "Factor2", "Factor3")], n_profiles=4)
lpa = estimate_profiles(policy_data[, c("Cigarette_tax", "Wage", "Parental_leave_type_cat", "Tax_cumulative_cat")], n_profiles=3)
lpa
plot_profiles(lpa)

policy_data$Class = get_data(lpa)$Class

boxplot(policy_data$IMR_2018 ~ policy_data$Class, main="A)", xlab="Latent Policy Profile", ylab="Infant Mortality Rate per 1,000")
policy_data$FMR_2014_18 = policy_data$Fetal_deaths_2018/policy_data$Births_2018*1000
boxplot(policy_data$FMR_2014_18 ~ policy_data$Class, main="B)", xlab="Latent Policy Profile", ylab="Fetal Mortality Rate per 1,000")

#export
write.csv(policy_data, file="data.csv", na="", row.names=F)


### CHOROPLETH US MAP ###

#retain only 50 U.S. states
us_carto = us_carto[us_carto$NAME %in% state.name, ]

#set projected coordinate system for U.S.
us_carto_proj = spTransform(us_carto,CRS("+init=epsg:2163"))

#need to transform AK, HI to fit under U.S. for single map; see https://stackoverflow.com/questions/13757771/relocating-alaska-and-hawaii-on-thematic-map-of-the-usa-with-ggplot2
fixup <- function(usa,alaskaFix,hawaiiFix){
  
  alaska=usa[usa$NAME=="Alaska",]
  alaska = fix1(alaska,alaskaFix)
  proj4string(alaska) <- proj4string(usa)
  
  hawaii = usa[usa$NAME=="Hawaii",]
  hawaii = fix1(hawaii,hawaiiFix)
  proj4string(hawaii) <- proj4string(usa)
  
  usa = usa[! usa$NAME %in% c("Alaska","Hawaii"),]
  usa = rbind(usa,alaska,hawaii)
  
  return(usa)
  
}

fix1 <- function(object,params){
  r=params[1];scale=params[2];shift=params[3:4]
  object = elide(object,rotate=r)
  size = max(apply(bbox(object),1,diff))/scale
  object = elide(object,scale=size)
  object = elide(object,shift=shift)
  object
}

us_map = fixup(us_carto_proj,c(-35,2,-2500000,-2500000),c(-35,1,5500000,-1600000))
rm(fix1,fixup,us_carto,us_carto_proj)

#choropleth shading for IMR
map_col = grey(seq(0.4, 0.9, by=0.25))
map_col_index = policy_data$Class[match(us_map$NAME, policy_data$State)]

#draw choropleth map
par(mar=rep(0.1,4))
plot(us_map,col=map_col[map_col_index])

legend("bottomright", legend=c("#1", "#2", "#3"), title="Latent Policy Profile", fill=map_col, cex=0.8, horiz=T)


### DESCRIPTIVES ###

describe(policy_data$Cigarette_tax[policy_data$Abbreviation!="DC"])
CrossTable(policy_data$Wage_fed_comparison_7.25[policy_data$Abbreviation!="DC"])
describe(policy_data$Wage[policy_data$Wage_fed_comparison_7.25=="Above" & policy_data$Abbreviation!="DC"])
CrossTable(policy_data$Parental_leave_type[policy_data$Abbreviation!="DC"])
CrossTable(policy_data$Tax_EITC[policy_data$Abbreviation!="DC"])
CrossTable(policy_data$Tax_CTC[policy_data$Abbreviation!="DC"])
CrossTable(policy_data$Tax_CDCTC[policy_data$Abbreviation!="DC"])
CrossTable(policy_data$Tax_cumulative_cat[policy_data$Abbreviation!="DC"])
