all: outline manuscript

outline: notes/outline.docx

manuscript: ms/manuscript.docx

notes/outline.docx: notes/outline.md template/reference.docx
	pandoc --reference-doc=template/reference.docx \
		--from=markdown \
		--to=docx \
		--csl=template/genome-research.csl \
		--bibliography=bib/references.bib \
		--output=notes/outline.docx \
		notes/outline.md

ms/manuscript.docx: ms/manuscript.md template/reference.docx
	pandoc --reference-doc=template/reference.docx \
		--from=markdown \
		--to=docx \
		--csl=template/genome-research.csl \
		--bibliography=bib/references.bib \
		--output=ms/manuscript.docx \
		ms/manuscript.md
