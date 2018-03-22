all: outline manuscript

outline: docx/outline.docx

manuscript: docx/manuscript.docx

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

#fig/fig-pc5.svg: src/fig-pc5.R
#	cd src; R CMD BATCH fig-pc5.R

fig/%.svg: src/%.R
	cd $(<D); Rscript --vanilla $(<F)
