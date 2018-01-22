all: outline

outline: notes/outline.docx

notes/outline.docx: notes/outline.md template/reference.docx
	pandoc --reference-docx=template/reference.docx \
		--from=markdown \
		--to=docx \
		--output=notes/outline.docx \
		notes/outline.md
