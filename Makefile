all: outline manuscript

outline: docx/outline.docx

manuscript: docx/manuscript.docx html/manuscript.html

docx/outline.docx: notes/outline.md template/reference.docx
	pandoc --reference-doc=template/reference.docx \
		--from=markdown \
		--to=docx \
		--csl=template/genome-research.csl \
		--bibliography=bib/references.bib \
		--output=docx/outline.docx \
		notes/outline.md

docx/manuscript.docx: ms/manuscript.md template/reference.docx fig/fig-pc5.svg
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
# css from https://gist.github.com/killercup/5917178
html/manuscript.html: ms/manuscript.md template/reference.docx fig/fig-pc5.svg tables/table-pc5-mapman.svg
	pandoc --from=markdown \
		--to=html \
		--css=../css/pandoc.css \
		--csl=template/genome-research.csl \
		--bibliography=bib/references.bib \
		--output=html/manuscript.html \
		ms/manuscript.md

fig/%.svg: src/%.R data/pca-rlog.Rdata
	cd $(<D); Rscript --vanilla $(<F)

tables/%.svg: src/%.R data/pca-rlog.Rdata data/mapman.Rdata
	cd $(<D); Rscript --vanilla $(<F)

data/pca-rlog.Rdata: src/get-rlog-pca.R data/dds.Rds
	cd src; Rscript --vanilla get-rlog-pca.R

data/mapman.Rdata: src/get-mapman.R data/MAPMAN\ BIN-Osa_MSU_v7.xlsx
