#### load ####
library(phyloseq)
library(dplyr)
library(tidyr)
library(qiime2R)
library(vegan)
library(readr)
library(tibble)

# Configurar el directorio de trabajo
setwd("~/Library/CloudStorage/GoogleDrive-danielavargasrobles@gmail.com/My Drive/Filipa_Daniela_2021 -present/Proyectos/VAGINAL msystems/ITS")
# Cargar datos de ASV
phy <- qza_to_phyloseq(features = "QIIME2/dada2_paired_end_table.qza")
ASV <- otu_table(phy)
data<-as.data.frame(t(ASV))
dim(data)

# Cargar datos de taxonomía
#new seqs
taxonomy <- read_tsv('QIIME2/exported-taxonomy/taxonomy.tsv', show_col_types = FALSE)

taxonomy$Taxon[1]

#gg_ext tambien funciona para eHOMD
taxonomy_cleaned0 <- taxonomy %>%
  dplyr::select(-Confidence) %>% # Eliminar la columna Confidence
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","SH"), sep = ";", fill = "right") %>%
  mutate(across(everything(), ~ gsub("^\\s+", "", .))) %>%
  mutate(across(everything(), ~ gsub(" ", "_", .))) %>%  # Reemplazar espacios por guiones bajos
  mutate(
    Kingdom = ifelse(is.na(Kingdom) | Kingdom == '' | Kingdom == 'k__', 'k__Unclassified', Kingdom),
    Phylum = ifelse(is.na(Phylum) | Phylum == '' | Phylum == 'p__', 'p__Unclassified', Phylum),
    Class = ifelse(is.na(Class) | Class == '' | Class == 'c__', 'c__Unclassified', Class),
    Order = ifelse(is.na(Order) | Order == '' | Order == 'o__', 'o__Unclassified', Order),
    Family = ifelse(is.na(Family) | Family == '' | Family == 'f__', 'f__Unclassified', Family),
    Genus = ifelse(is.na(Genus) | Genus == '' | Genus == 'g__', 'g__Unclassified', Genus),
    Species = ifelse(is.na(Species) | Species == '' | Species == 's__', 's__Unclassified', Species),
    Species = ifelse(Genus != 'g__Unclassified' & Species == 's__Unclassified', 's__Unclassified', 
                     ifelse(Species != 's__Unclassified', sub('s__', '', Species), 
                            "g__Unclassified_s__Unclassified"))
  ) %>%
  column_to_rownames(var = "Feature ID") %>%  # Establecer Feature ID como nombres de fila
  filter(Kingdom != "k__Archaea" , Kingdom!= "k__Eukaryota_kgd_Incertae_sedis",Kingdom != "k__Protista", Kingdom != "Unassigned" , Phylum != "p__Unclassified", 
         Class != "c__Unclassified", Order != "o__Unclassified", Family!= "f__mitochondria",)

# Supongamos que 'taxonomy_cleaned0' es un vector o data frame que contiene las taxonomías



head(taxonomy_cleaned0)

taxonomy_cleaned <- tax_table(as.matrix(taxonomy_cleaned0))

# PRUNE taxonomy after ASV table filtering
common_taxa <- intersect(taxa_names(taxonomy_cleaned), taxa_names(ASV))
TAXA_clean <- prune_taxa(common_taxa, taxonomy_cleaned)
OTU_clean <- prune_taxa(common_taxa, ASV)

# Create the phyloseq object
phy1 <- phyloseq(TAXA_clean, OTU_clean)
phy1


#### Filtrar taxa por abundancia total > 1
phy_filtered <- prune_taxa(taxa_sums(phy1) > 0, phy1) 

# Calcular la presencia (número de muestras en las que cada ASV está presente)
asv_presence <- apply(otu_table(phy_filtered), 1, function(x) sum(x > 0))

# Selectiona los ASVs que están presentes en 2 o mas muestras
phy1 <- prune_taxa(asv_presence > 0, phy_filtered)
phy1

# Extraer la tabla de taxonomía del objeto phyloseq
tax_table <- as.data.frame(tax_table(phy1))



#####para verificar qué no está clasificado

# Filtrar las ASVs que no están clasificadas a nivel de especie
unclassified_asvs <- tax_table(phy1) %>%
  as.data.frame() %>%
  #filter(grepl("Unclassified", Species)) %>%
  #filter(grepl("Fungi_sp", Species)) %>%
  filter(grepl("Unclassified", Genus)) %>%
  
  rownames()

# Extraer la tabla de OTU del objeto phyloseq
otu_table <- otu_table(phy1)

# Filtrar la tabla de OTU para quedarse solo con las ASVs sin clasificar en especie
unclassified_otu_table <- otu_table[unclassified_asvs, ]

# Extraer la información taxonómica (familia y género)
unclassified_taxonomy <- tax_table(phy1) %>%
  as.data.frame() %>%
  filter(rownames(.) %in% unclassified_asvs) %>%
  dplyr::select(Family, Genus)

# Calcular el número total de lecturas (reads) para cada ASV sin clasificar
unclassified_counts <- rowSums(unclassified_otu_table)

# Calcular el número de muestras que contienen cada ASV sin clasificar
unclassified_sample_counts <- apply(unclassified_otu_table, 1, function(x) sum(x > 0))

# Crear un dataframe con los resultados, incluyendo la familia y el género
unclassified_summary <- data.frame(
  ASV_ID = unclassified_asvs,
  Family = unclassified_taxonomy$Family,
  Genus = unclassified_taxonomy$Genus,
  Total_Reads = unclassified_counts,
  Sample_Count = unclassified_sample_counts
)

# Ordenar los resultados por el número total de lecturas (descendente)
unclassified_summary %>%
  arrange(desc(Total_Reads))


# # Lista de ASV_IDs específicos
 ASV_IDs <- c(
   "e798ca7083a01a4e8107375f1a3f3891",
   "fc4fd617c411600ade080e8e3da15c59",
   "ec0ad58cb5a75324cea77a6a92225664",
   "08593fe3dba5bd28fe325c0e71297e73",
   "4f4802ffc515544754e84d1c36c9bb00",
   "a4e416936faa6f84dfca6d10059df489",
   "e30b933396c0faf703f2d6ef46b63db0",
   "ec8921fa949341e2c955588d758c9ec0",
   "0544621590e254d7c4f7089f036ed640",
   "32ee9be0fd8de61b65377d2b827d977c",
   "bd17952012d583987631b7ba47a7a106",
   'b216dc31242e2e3ce3681758474e049b',
   "74b880da69d82efab711e88890ecd974",
   "698a82494a97a929432405646aeabff9",
   "55f794b5e841c6e1e3f49e4a139494b7",
   "02ae2ddc274f1fa06694b08a4c4abe95",
   "58e0660790d5cb51697bede87ce67f69",
   "882bee54a0b623f1f02ebc6285476903",
   "197ba64bbe986021020d6ad9e55881f0",
   "1ac1428979175bf2857e8b785eb4f746",
   "570854ea9fabdf76148ffcdb90d09da2",
   "75666fdcb60c738f01a2de8494aafd6f",
   "cb11edb49620f485aa3956c64764dcb3"
 )
 
# # Extraer la tabla de taxonomía para los ASV específicos
 taxonomy_subset <- tax_table(phy1)[ASV_IDs, ]
 taxonomy_subset
 write.csv(taxonomy_subset, "g__Unclassified", row.names = TRUE)
 
 # Lista de nombres de ASV de interés g__Fungi_gen_Incertae_sedis
 ASV_IDs <- c("f7c6544988fd56fdc6614d8b67bc0ecd", "eccaadc7974cbe1f434a3017f470f90c", 
               "004589b84321c4f78fe4d2164ff8e567", "17994b3a7615f5112d10f30c5ac26669", 
               "0bd921e085b3e7b1377f5c2e597ee023", "79512e8a40a58b1bd9c8dbbc7a9167b6", 
               "5f158cac8629072384c6402f8cf7b8f4", "0bdc52ff2a1f11d597ba1c6f9d1b4ef9", 
               "909b8a4885eb8cb90465376d54a0fd67", "66707192a63e13b02bbf997cf5234084", 
               "d6618c2440665f462e13a0d0b81379d7", "73aafae4fb1166c87b2846a8836b7a00", 
               "691304e567f004b35e2a54b6353b3c64", "d5981d4942d5fb8e39719e56eec16277", 
               "0b6c910f35e1aeaa4b3a71c61c89ef42", "b3978ade5084e2d48f826cdc25db4c4d", 
               "224ef90e5cf34059addab2df2a67f22d", "3b982659f2b4c91a5208d04801425896", 
               "ca52c2120a4cf5e2d0bac59ac48edddc", "156bb949755dd52447dcae80a18f537a", 
               "bcc700a97077b12517c42adb4b28af2a", "8e878d4c34b2e8a0cd10840c752da943", 
               "83551a3808bb06eee684702d96647d8b", "c80ba35892338c644c36aead29467d59", 
               "22144b1aa70a59ecaa326cd4e147be22", "59c83662ca40b61324a0b645231938d2", 
               "1d1a0b8fb8ae9dbf17b8d7cc03c1e9b3", "7cef19db3c8f4cabcc4d9c4ab0ea5328", 
               "f0dc87d13ceb6230853d4cc9f4098c28", "e23dda43b76ebbf4a16ef45325673e14", 
               "8bc194d488378d7c4cb1739579ff6067", "8ddba597f29d43b78e476e40e6971536", 
               "5272b6c55febcc1b737eb54f01e375cf", "5c9ece9b1928e1ade960fbe6d9b9b186", 
               "e2556a2b20c7a3a059ae543dc2007aac", "d5562305ce756f0f6760bb2370444383", 
               "4aeb882d1bbf8388ac02c4d89c634178", "bb6edf13a734542db4e6a22859bcea1c", 
               "9703e8eba18cf52bd9656389738bfc80", "d189d21ae4da43736e71534095a3a320", 
               "025d570cddc3562472f2e68f28c2c9ef", "6de4c33257a9b51646ea3a19471cb975", 
               "9097dad2f00593d16c34c69e52145e40", "8252fb76f79ab2d1991c8a5619fbe728", 
               "ac9d4910a27f0b12b3bf9ac47b8e44fe")
 
 # # Extraer la tabla de taxonomía para los ASV específicos
 taxonomy_subset <- tax_table(phy1)[ASV_IDs, ]
 taxonomy_subset
 write.csv(taxonomy_subset, "Fungi_gen_Incertae_sedis", row.names = TRUE)
 
 
 
 
 
#################### METADATA ###########################

 # Cargar la metadata
 met = read.csv("/Users/danielavargasrobles/Library/CloudStorage/GoogleDrive-danielavargasrobles@gmail.com/My Drive/Filipa_Daniela_2021 -present/Proyectos/VAGINAL msystems/REANALYSIS/met/clean_errata_walpha.csv")
 
 # Corregir la columna Lesion
 met <- met %>%
   mutate(Lesion = case_when(
     is.na(Lesion) & consensus.worked_no_ascus_chSCCtoHGSIL == "HGSIL" ~ "Positive",
     is.na(Lesion) & consensus.worked_no_ascus_chSCCtoHGSIL == "LGSIL" ~ "Positive",
     is.na(Lesion) & consensus.worked_no_ascus_chSCCtoHGSIL == "NILM" ~ "Negative",
     Lesion == "Negative" & consensus.worked_no_ascus_chSCCtoHGSIL == "HGSIL" ~ "Positive",
     Lesion == "Negative" & consensus.worked_no_ascus_chSCCtoHGSIL == "LGSIL" ~ "Positive",
     Lesion == "Positive" & consensus.worked_no_ascus_chSCCtoHGSIL == "NILM" ~ "Negative",
     Lesion == "Negative" & consensus.worked_no_ascus_chSCCtoHGSIL == "NILM" ~ "Negative",
     TRUE ~ Lesion  # Retains the original value if it is not "Positive"
   ))
 
 # Mostrar las filas seleccionadas
 met %>%
   dplyr::select(Lesion, consensus.worked_no_ascus_chSCCtoHGSIL, SubjectID) %>%
   arrange(Lesion)
 
 
 
table(met$Lesion, met$consensus.worked_no_ascus_chSCCtoHGSIL)
TAXA_clean=as.matrix(tax_table(phy1))
OTU_clean=as.matrix(otu_table(phy1, taxa_are_rows = TRUE))

# Remover "_ITS" o "ITS_" de los nombres de las columnas
colnames(OTU_clean) <- gsub("_ITS|ITS_", "", colnames(OTU_clean))
phy1 <- phyloseq(otu_table(OTU_clean, taxa_are_rows = TRUE), 
                 tax_table(TAXA_clean))

# Visualizar el resultado
head(OTU_clean)
dim(OTU_clean)

# Reordenar metadata para que coincida con la tabla ASV
met_ordered <- met[match(colnames(OTU_clean), met$SampleID), ]
if(any(is.na(met_ordered$SampleID))){
  stop("Some samples in OTU_clean do not match with metadata")
}


# Asignar los nombres de las columnas de filter_otu a las filas de met_ordered
rownames(met_ordered) <- colnames(OTU_clean)

#check
head(colnames(OTU_clean))
head(met_ordered$SampleID)

# Añadir la metadata al objeto phyloseq
rownames(met_ordered) <- sample_names(phy1)
met_ordered1 <- sample_data(met_ordered)
dim(TAXA_clean)

PHYraw <- phyloseq(TAXA_clean, OTU_clean, met_ordered1)

PHYraw

save(PHYraw,file="R/phy_objects/PHYraw.RData")
load("R/phy_objects/PHYraw.RData")


##### Filtering #####

# Calcular la suma de secuencias por muestra
sumatoria <- sample_sums(PHYraw)
min(sumatoria)

sumatoria_ordenada <- sort(sumatoria)
print(head(sumatoria_ordenada, 30))

# Agregar la suma de secuencias como una nueva columna en la metadata
sample_data(PHYraw)$num_seqs <- sumatoria

seq_to_check=data.frame(sample_data(PHYraw))


seq_to_check %>% 
  dplyr::select(num_seqs, SubjectID, Visit_number) %>% 
  arrange(Visit_number,-num_seqs)


# Filtrar las muestras con menos de 100 seqs
PHYraw

seq_to_check=data.frame(sample_data(PHYraw))
seq_to_check %>% 
  dplyr::select(num_seqs, SubjectID, Visit_number) %>% 
  arrange(SubjectID,Visit_number,-num_seqs)

# Crear una lista de SampleID para eliminar porque estas muestrsa están tienen dos timepoints. De los duplicados se está eliminanado el 2_visit 
samples_to_remove <- c("2156b", "2193b", "2206b", "2217b", "2274b") #esta bien, despues de remover las muestrs con menos de  100 seqs
PHYraw <- subset_samples(PHYraw, !(SampleID %in% samples_to_remove))
PHYraw

#Filtrar muestras con menos de 100 secuencias
PHYraw <- prune_samples(sample_sums(PHYraw) > 100, PHYraw)
#eliminar de phy las muestras de shannon_ITS == NA
PHYraw <- prune_samples(sample_sums(PHYraw) > 0, PHYraw)
# Eliminar taxones sin presencia en el objeto phy filtrado
PHYraw <- prune_taxa(taxa_sums(PHYraw) > 0, PHYraw)
PHYraw

#### Phyloseq at ASV####
# 1. Extraer tabla ASV 
OTU1 = as(otu_table(PHYraw), "matrix")
OTUdf = as.data.frame(OTU1)

# Calculando las 4 métricas de diversidad alfa a nivel ASV
# Calculando Shannon y Simpson, pero con una condición para cuando obs_ITS es 0

obs = apply(OTUdf, 2, function(x){
  nws = length(which(x > 0))
  return(nws)
})

# Reemplazamos los cálculos de Shannon y Simpson con condiciones
shannon = apply(OTUdf, 2, function(x) {
  if (sum(x) == 0) {
    return(NA)  # Si no hay lecturas, asignar NA
  } else {
    return(diversity(x, index = "shannon", base = exp(1)))  # Calcular Shannon solo si hay lecturas
  }
})

simpson = apply(OTUdf, 2, function(x) {
  if (sum(x) == 0) {
    return(NA)  # Si no hay lecturas, asignar NA
  } else {
    return(diversity(x, index = "simpson"))  # Calcular Simpson solo si hay lecturas
  }
})

# Calcular Chao1 usando estimate_richness
chao = estimate_richness(PHYraw, split = TRUE, measures = "Chao1")

# Crear dataframe con los resultados
div_asv = data.frame(shannon_ITS = shannon, 
                     simpson_ITS = simpson,
                    # chao1_ITS = chao,
                     obs_ITS = obs)

# Reemplazar Shannon y Simpson con NA cuando el número de lecturas es 0
div_asv$shannon_ITS[div_asv$obs_ITS == 0] <- NA
div_asv$simpson_ITS[div_asv$obs_ITS == 0] <- NA

#merging columns of alpha metrics with the rest of the metadata by SampleID
md1 = sample_data(PHYraw)

asv=cbind(md1,div_asv)
rownames(asv) = sample_names(PHYraw)
sampledata = sample_data(asv)

#Re doing phyloseq object with new added columns (alpha div metrics) in metadata
PHYrawITS_alpha = phyloseq(otu_table(PHYraw),tax_table(PHYraw), sampledata)#phy_tree(PHYraw) )
hola=sample_data(PHYrawITS_alpha)
save(PHYrawITS_alpha,file="R/phy_objects/PHYrawITS_alpha.RData")
load("R/phy_objects/PHYrawITS_alpha.RData")


####Phyloseq at genus ####
# 1. Agrupar a nivel de género
PHYgenus <- tax_glom(PHYraw, taxrank = "Genus")

# Extraer la tabla de ASV agrupada a nivel de género
OTU_genus = as(otu_table(PHYgenus), "matrix")
OTUdf_genus = as.data.frame(OTU_genus)

# Calculando las métricas de diversidad alfa a nivel de género
# Calculando Shannon y Simpson, con una condición para cuando obs_ITS es 0

obs_genus = apply(OTUdf_genus, 2, function(x){
  nws = length(which(x > 0))
  return(nws)
})

# Reemplazamos los cálculos de Shannon y Simpson con condiciones
shannon_genus = apply(OTUdf_genus, 2, function(x) {
  if (sum(x) == 0) {
    return(NA)  # Si no hay lecturas, asignar NA
  } else {
    return(diversity(x, index = "shannon", base = exp(1)))  # Calcular Shannon solo si hay lecturas
  }
})

simpson_genus = apply(OTUdf_genus, 2, function(x) {
  if (sum(x) == 0) {
    return(NA)  # Si no hay lecturas, asignar NA
  } else {
    return(diversity(x, index = "simpson"))  # Calcular Simpson solo si hay lecturas
  }
})

# Calcular Chao1 a nivel de género usando estimate_richness
chao_genus = estimate_richness(PHYgenus, split = TRUE, measures = "Chao1")

# Crear dataframe con los resultados de diversidad alfa a nivel de género
div_genus = data.frame(shannon_genus = shannon_genus, 
                       simpson_genus = simpson_genus,
                    #   chao1_genus = chao_genus,
                       obs_genus = obs_genus)

# Reemplazar Shannon y Simpson con NA cuando el número de lecturas es 0
div_genus$shannon_genus[div_genus$obs_genus == 0] <- NA
div_genus$simpson_genus[div_genus$obs_genus == 0] <- NA

# Integrar las métricas de diversidad alfa con la metadata por SampleID
md_genus = sample_data(PHYgenus)
genus_alpha_data = cbind(md_genus, div_genus)
rownames(genus_alpha_data) = sample_names(PHYgenus)
sampledata_genus = sample_data(genus_alpha_data)


# Rehacer el objeto phyloseq con las nuevas columnas de diversidad alfa en la metadata
PHYraw_genus_alpha_ITS = phyloseq(otu_table(PHYgenus), tax_table(PHYgenus), sampledata_genus)

#extrar metedatad as dataframe
met=as.data.frame(sample_data(PHYraw_genus_alpha_ITS))

table(met$consensus.worked_no_ascus_chSCCtoHGSIL)
# Guardar el objeto actualizado con métricas de diversidad alfa a nivel de género
save(PHYraw_genus_alpha_ITS, file = "R/phy_objects/PHYraw_genus_alpha_ITS.RData")
load("R/phy_objects/PHYraw_genus_alpha_ITS.RData")


#### IMPUTAR 1 BMI ####
# Convertir los metadatos del objeto phyloseq en un data.frame
met <- data.frame(sample_data(PHYraw_genus_alpha_ITS))

met %>% 
  dplyr::select(SampleID, Group, BMI) %>% 
  arrange( BMI)
dim(met)
#al final no habia muestra para imputar
# if ("BMI" %in% colnames(met) && any(is.na(met$BMI))) {
#   met$BMI[is.na(met$BMI)] <- tapply(met$BMI, met$Group, median, na.rm = TRUE)[met$Group[is.na(met$BMI)]]
#   sample_data(PHYraw_genus_alpha_ITS) <- sample_data(met)
# }

PHYraw_genus_alpha_ITS <- prune_samples(sample_sums(PHYraw_genus_alpha_ITS) > 0, PHYraw_genus_alpha_ITS)
PHYraw_genus_alpha_ITS <- prune_taxa(taxa_sums(PHYraw_genus_alpha_ITS) > 0, PHYraw_genus_alpha_ITS)
PHYraw_genus_alpha_ITS

save(PHYraw_genus_alpha_ITS, file = "R/phy_objects/PHYraw_genus_alpha_ITS.RData")

### Save PHY at genus ####

load("R/phy_objects/PHYraw_genus_alpha_ITS.RData")




# load PHYrawITS_alpha para hacer el mismo filtrado
load("R/phy_objects/PHYrawITS_alpha.RData")
PHYrawITS_alpha <- subset_samples(PHYrawITS_alpha, !(SampleID %in% samples_to_remove))
PHYrawITS_alpha <- prune_samples(sample_sums(PHYrawITS_alpha) > 0, PHYrawITS_alpha)
PHYrawITS_alpha <- prune_taxa(taxa_sums(PHYrawITS_alpha) > 0, PHYrawITS_alpha)
PHYrawITS_alpha

#IMPUTAR 1 BMI
# Convertir los metadatos del objeto phyloseq en un data.frame
met <- data.frame(sample_data(PHYrawITS_alpha))
# 
# if ("BMI" %in% colnames(met) && any(is.na(met$BMI))) {
#   met$BMI[is.na(met$BMI)] <- tapply(met$BMI, met$Group, median, na.rm = TRUE)[met$Group[is.na(met$BMI)]]
#   sample_data(PHYrawITS_alpha) <- sample_data(met)
# }

save(PHYrawITS_alpha, file = "R/phy_objects/PHYrawITS_alpha.RData")

load("R/phy_objects/PHYrawITS_alpha.RData")

#extrat metadata as dataframe
met=as.data.frame(sample_data(PHYrawITS_alpha))
table(met$consensus.worked_no_ascus_chSCCtoHGSIL)
############################### PHY 12 sp ######################
species_to_keep <- c("Candida_albicans", "Candida_tropicalis", "Candida_nivariensis",
                     "Candida_parapsilosis", "Hortaea_werneckii", "Saccharomyces_cerevisiae",
                     "Nakaseomyces_sp", "Candida_metapsilosis", "Trichosporon_asahii",
                     "Meyerozyma_caribbica", "Cystobasidium_calyptogenae", "Malassezia_restricta")

# Filtra el objeto phyloseq
phy_filtered <- subset_taxa(PHYrawITS_alpha, Species %in% species_to_keep)

#verificar las Species que quedaron la columna Species del taxa table
tax_table(phy_filtered)
phy_filtered

#extrat metadata as dataframe
met=as.data.frame(sample_data(phy_filtered))
table(met$consensus.worked_no_ascus_chSCCtoHGSIL)
#COmo hay tantos ASV para cada species vamos a aglomerar
# Aglomerar los ASV a nivel de especie
phy <- tax_glom(phy_filtered, taxrank = "Species")
phy


met=as.data.frame(sample_data(phy))
table(met$consensus.worked_no_ascus_chSCCtoHGSIL)


#Calcular nuvos indices dde shannon obs y simpson y meterlos en ela metadata
alpha_diversity <- estimate_richness(phy, measures = c("Observed", "Shannon", "Simpson"))
#guardar estos indicas como columas en la metadata del objeto phyloseq
sample_data(phy)$obs_ITS_12sp <- alpha_diversity$Observed
sample_data(phy)$shannon_ITS_12sp <- alpha_diversity$Shannon
sample_data(phy)$simpson_ITS_12sp <- alpha_diversity$Simpson

# remove samples que se hayan quedado sin taxones
PHYrawITS_alpha_12sp <- prune_samples(sample_sums(phy) > 0, phy)
PHYrawITS_alpha_12sp



# Guardar el objeto phyloseq filtrado R/phy_objects/PHYrawITS_alpha_12sp.RData
save(PHYrawITS_alpha_12sp, file = "R/phy_objects/PHYrawITS_alpha_12sp.RData")





#### filtering ASV in 2 samples####

load("R/phy_objects/PHYraw.RData")
PHYraw1<-PHYraw

# Filtrar las muestras con menos de 100 seqs
PHYraw1 <- prune_samples(sample_sums(PHYraw1) > 100, PHYraw1)
seq_to_check=data.frame(sample_data(PHYraw1))

# Crear una lista de SampleID para eliminar porque estas muestrsa están tienen dos timepoints. De los duplicados se está eliminanado el 2_visit 
samples_to_remove <- c("2156b", "2193b", "2206b", "2217b", "2274b") #esta bien, despues de remover las muestrs con menos de  100 seqs
PHYraw1 <- subset_samples(PHYraw1, !(SampleID %in% samples_to_remove))
PHYraw1
#eliminar de phy las muestras de shannon_ITS == NA
PHYraw1 <- prune_samples(sample_sums(PHYraw1) > 0, PHYraw1)
PHYraw1
#save
save(PHYraw1,file="R/phy_objects/PHYraw1.RData")
PHYraw2 <- prune_taxa(rowSums(otu_table(PHYraw1) > 0) >= 2, PHYraw1)
PHYraw2

##ALPHA diversity indices##
# 1. Extraer tabla ASV 
OTU1 = as(otu_table(PHYraw2), "matrix")
OTUdf = as.data.frame(OTU1)

# Calculando las 4 métricas de diversidad alfa a nivel ASV
# Calculando Shannon y Simpson, pero con una condición para cuando obs_ITS es 0

obs = apply(OTUdf, 2, function(x){
  nws = length(which(x > 0))
  return(nws)
})

# Reemplazamos los cálculos de Shannon y Simpson con condiciones
shannon = apply(OTUdf, 2, function(x) {
  if (sum(x) == 0) {
    return(NA)  # Si no hay lecturas, asignar NA
  } else {
    return(diversity(x, index = "shannon", base = exp(1)))  # Calcular Shannon solo si hay lecturas
  }
})

simpson = apply(OTUdf, 2, function(x) {
  if (sum(x) == 0) {
    return(NA)  # Si no hay lecturas, asignar NA
  } else {
    return(diversity(x, index = "simpson"))  # Calcular Simpson solo si hay lecturas
  }
})

# Calcular Chao1 usando estimate_richness
chao = estimate_richness(PHYraw2, split = TRUE, measures = "Chao1")

# Crear dataframe con los resultados
div_asv = data.frame(shannon_ITS = shannon, 
                     simpson_ITS = simpson,
                     # chao1_ITS = chao,
                     obs_ITS = obs)

# Reemplazar Shannon y Simpson con NA cuando el número de lecturas es 0
div_asv$shannon_ITS[div_asv$obs_ITS == 0] <- NA
div_asv$simpson_ITS[div_asv$obs_ITS == 0] <- NA

#merging columns of alpha metrics with the rest of the metadata by SampleID
md1 = sample_data(PHYraw2)

asv=cbind(md1,div_asv)
rownames(asv) = sample_names(PHYraw2)
sampledata = sample_data(asv)


#Re doing phyloseq object with new added columns (alpha div metrics) in metadata
PHYraw_2_ITS_alpha = phyloseq(otu_table(PHYraw2),tax_table(PHYraw2), sampledata)#phy_tree(PHYraw) )
hola=sample_data(PHYraw_2_ITS_alpha)
PHYraw_2_ITS_alpha
save(PHYraw_2_ITS_alpha,file="R/phy_objects/PHYraw_2_ITS_alpha.RData")



### SH level####
#usamos PHYraw1.RData porque solo quiero eliminar las muestras menores de 100 secuencias pero nNO eliminar los ASV que no estan presentes en menos de 2 muestras

load("R/phy_objects/PHYraw1.RData")
PHYrawSH<-PHYraw1
#aglomerar a SH level
taxa <- data.frame(tax_table(PHYrawSH))
taxa$

PHYrawSH <- tax_glom(PHYrawSH, taxrank = "SH")
#remove SH with less than 2 sequences
PHYrawSH <- prune_taxa(taxa_sums(PHYrawSH) > 1, PHYrawSH)
taxa <- data.frame(tax_table(PHYrawSH))
taxa
##ALPHA diversity indices##
# 1. Extraer tabla ASV 
OTU1 = as(otu_table(PHYrawSH), "matrix")
OTUdf = as.data.frame(OTU1)

# Calculando las 4 métricas de diversidad alfa a nivel ASV
# Calculando Shannon y Simpson, pero con una condición para cuando obs_ITS es 0

obs = apply(OTUdf, 2, function(x){
  nws = length(which(x > 0))
  return(nws)
})

# Reemplazamos los cálculos de Shannon y Simpson con condiciones
shannon = apply(OTUdf, 2, function(x) {
  if (sum(x) == 0) {
    return(NA)  # Si no hay lecturas, asignar NA
  } else {
    return(diversity(x, index = "shannon", base = exp(1)))  # Calcular Shannon solo si hay lecturas
  }
})

simpson = apply(OTUdf, 2, function(x) {
  if (sum(x) == 0) {
    return(NA)  # Si no hay lecturas, asignar NA
  } else {
    return(diversity(x, index = "simpson"))  # Calcular Simpson solo si hay lecturas
  }
})

# Calcular Chao1 usando estimate_richness
chao = estimate_richness(PHYrawSH, split = TRUE, measures = "Chao1")

# Crear dataframe con los resultados
div_asv = data.frame(shannon_ITS = shannon, 
                     simpson_ITS = simpson,
                     # chao1_ITS = chao,
                     obs_ITS = obs)

# Reemplazar Shannon y Simpson con NA cuando el número de lecturas es 0
div_asv$shannon_ITS[div_asv$obs_ITS == 0] <- NA
div_asv$simpson_ITS[div_asv$obs_ITS == 0] <- NA

#merging columns of alpha metrics with the rest of the metadata by SampleID
md1 = sample_data(PHYrawSH)

asv=cbind(md1,div_asv)
rownames(asv) = sample_names(PHYrawSH)
sampledata = sample_data(asv)


#Re doing phyloseq object with new added columns (alpha div metrics) in metadata
PHYraw_SH_ITS_alpha = phyloseq(otu_table(PHYrawSH),tax_table(PHYrawSH), sampledata)#phy_tree(PHYraw) )
hola=sample_data(PHYraw_SH_ITS_alpha)
PHYraw_SH_ITS_alpha
save(PHYraw_SH_ITS_alpha,file="R/phy_objects/PHYraw_SH_ITS_alpha.RData")




# 2 muestras sin SH ####
#explorando qué clasificacion taxonomica tienen las dos muestras que no tienen classificacion en SH

phy_2128b <- prune_samples(sample_names(PHYraw1) == "2128b", PHYraw1)

# Eliminar taxa con abundancia 0 en la muestra
phy_2128b_filtered <- prune_taxa(taxa_sums(phy_2128b) > 0, phy_2128b)

# Extraer la tabla taxonómica de la muestra filtrada
taxa_2128b <- as.data.frame(tax_table(phy_2128b_filtered)) %>%
  rownames_to_column("Feature_ID")  # Mantener identificador de OTUs/ASVs/SH

#ver del dataframe
taxa_2128b %>% 
  dplyr::select(Genus, Species, SH)
dim(taxa_2128b)


###ahora la otra: 2292

phy_2292 <- prune_samples(sample_names(PHYraw1) == "2292", PHYraw1)

# Eliminar taxa con abundancia 0 en la muestra
phy_2292_filtered <- prune_taxa(taxa_sums(phy_2292) > 0, phy_2292)

# Extraer la tabla taxonómica de la muestra filtrada
phy_2292 <- as.data.frame(tax_table(phy_2292_filtered)) %>%
  rownames_to_column("Feature_ID")  # Mantener identificador de OTUs/ASVs/SH

#ver del dataframe
phy_2292 %>% 
  dplyr::select(Family, Genus, SH)
dim(taxa_2128b)
