# The outline is just for us
all: outline manuscript

outline: docx/outline.docx

# The official manuscript is in docx, but the html version is much easier to use
manuscript: docx/manuscript.docx html/manuscript.html

# For simplicity, specify common dependencies for both the manuscript and the
# outline

TEXT_DEPEND = \
	template/reference.docx \
	fig/fig-branches-boxplot.svg \
	fig/fig-spn-prediction.svg \
	fig/fig-pc5.svg \
	tables/table-pc5-mapman.svg \
	bib/references.bib \
	template/genome-research.csl \
	css/pandoc.css

docx/outline.docx: notes/outline.md $(TEXT_DEPEND)
	pandoc --reference-doc=template/reference.docx \
		--from=markdown \
		--to=docx \
		--csl=template/genome-research.csl \
		--bibliography=bib/references.bib \
		--output=docx/outline.docx \
		notes/outline.md

docx/manuscript.docx: ms/manuscript.md $(TEXT_DEPEND)
	pandoc --reference-doc=template/reference.docx \
		--from=markdown \
		--to=docx \
		--csl=template/genome-research.csl \
		--bibliography=bib/references.bib \
		--output=docx/manuscript.docx \
		ms/manuscript.md

# Produce an html manuscript in the main folder
# the html manuscript is just for us
# It is much nicer to read than the docx
html/manuscript.html: ms/manuscript.md $(TEXT_DEPEND)
	pandoc --from=markdown \
		--to=html \
		--css=../css/pandoc.css \
		--csl=template/genome-research.csl \
		--bibliography=bib/references.bib \
		--output=html/manuscript.html \
		ms/manuscript.md

# For simplicity, specify common DATA dependecies for figures and tables
FIG_TAB_DEPEND = \
  data/phenotype-tidy.Rdata \
	data/pca-rlog.Rdata \
	data/mapman.Rdata \
	src/helper-functions.R

fig/%.svg: src/%.R $(FIG_TAB_DEPEND)
	cd $(<D); Rscript --vanilla $(<F)

tables/%.svg: src/%.R $(FIG_TAB_DEPEND)
	cd $(<D); Rscript --vanilla $(<F)

# Starting point of the analysis from which everything else depends

data/pca-rlog.Rdata: \
 src/get-rlog-pca.R data/dds.Rds
	cd $(<D); Rscript --vanilla $(<F)

data/mapman.Rdata: \
 src/get-mapman.R data/MAPMAN\ BIN-Osa_MSU_v7.xlsx \
 src/helper-functions.R
	cd $(<D); Rscript --vanilla $(<F)

data/phenotype-tidy.Rdata: \
 src/wrangle-phenotypes.R \
 data/Phenotype_PanicleSequenced_corrected.xlsx\
 data/OsOgObOrPTRAPdata_PaperTom.txt
	cd $(<D); Rscript --vanilla $(<F)
