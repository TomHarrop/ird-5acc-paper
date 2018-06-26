library(tidyverse)
library(biomaRt)


# Connect to gramene ------------------------------------------------------

# tried different options but I get curl errors

# listMarts(host = "gramene.org")
# gramene <- useMart(biomart = "ENSEMBL_MART_PLANT",
#                    host = "gramene.org")
# # listDatasets(gramene)
# 
# # dataset: japonica IRGSP-1.0
# # japonica <- useDataset(dataset = "osativa_eg_gene",
# #                        mart = gramene)
# sativa <- useMart(biomart = "ENSEMBL_MART_PLANT",
#                   host = "http://ensembl.gramene.org",
#                   dataset = "osativa_eg_gene")

# 


# View(listFilters(japonica))
# View(listAttributes(japonica))


# get MSU to RAPDB --------------------------------------------------------

dict_url <-  "http://rapdb.dna.affrc.go.jp/download/archive/RAP-MSU_2018-03-29.txt.gz"


dict_path <- "../data/msu-to-rapdb.Rdata"
if(!file.exists(dict_path)) {
  dict <- read_delim(dict_url,
                     delim = "\t", 
                     col_names = c("rapdb", "msu"))
  save(dict, file = dict_path)
} else {
  load(dict_path)
}  

id <- "LOC_Os03g12950"
id_rap <- dict %>%
  filter(grepl(id, msu)) %>%
  pull(rapdb)

# getBM(attributes = "ensembl_gene_id",
#       mart = sativa) 

# getBM(filters = "chromosome_name",
#       values = 1,
#       attributes = "ensembl_gene_id",
#       mart = sativa)


# Try direct XML query ----------------------------------------------------

tst <- getXML(host = "http://ensembl.gramene.org/biomart/martservice?",
              xmlquery = '<?xml version="1.0" encoding="UTF-8"?>
         <!DOCTYPE Query>
         <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
         
         <Dataset name = "osativa_eg_gene" interface = "default" >
         <Filter name = "chromosome_name" value = "1"/>
         <Attribute name = "ensembl_gene_id" />
         <Attribute name = "ensembl_transcript_id" />
         </Dataset>
         </Query>')

tst <- getXML(host = "http://ensembl.gramene.org/biomart/martservice?",
              xmlquery = '<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "osativa_eg_gene" interface = "default" >
		<Filter name = "upstream_flank" value = "2000"/>
		<Filter name = "ensembl_gene_id" value = "Os03g0232200"/>
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_transcript_id" />
		<Attribute name = "gene_flank" />
	</Dataset>
</Query>')
