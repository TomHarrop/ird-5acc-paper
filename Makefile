all: dot_points

dot_points: notes/dot_points.docx

notes/dot_points.docx: notes/dot_points.md template/reference.docx
	pandoc --reference-docx=template/reference.docx \
		--from=markdown \
		--to=docx \
		--output=notes/dot_points.docx \
		notes/dot_points.md
