manuscript: docx/manuscript.docx 
si: pdf/supporting_information.pdf docx/supporting_information.docx
ms_pdf: pdf/manuscript.pdf

docx/manuscript.docx: ms/front_matter.md ms/abstract.md ms/introduction.md ms/methods.md ms/results.md ms/discussion.md ms/end_matter.md template/ref_loc.md ms/figure_legends.md ms/si_list.md bib/references.bib template/reference.docx template/new-phytologist.csl template/ref_loc.md
	pandoc --reference-doc=template/reference.docx \
		--filter pandoc-citeproc \
		--from=markdown \
		--to=docx \
		--bibliography=bib/references.bib \
		--csl=template/new-phytologist.csl \
		-o docx/manuscript.docx \
		ms/front_matter.md \
		ms/abstract.md \
		ms/introduction.md \
		ms/methods.md \
		ms/results.md \
		ms/discussion.md \
		ms/end_matter.md \
		template/ref_loc.md \
		ms/figure_legends.md \
		ms/si_list.md

pdf/manuscript.pdf: ms/front_matter.md ms/abstract.md ms/introduction.md ms/methods.md ms/results.md ms/discussion.md ms/end_matter.md template/ref_loc.md ms/figure_legends.md ms/si_list.md bib/references.bib template/reference.docx template/new-phytologist.csl template/ref_loc.md
	pandoc \
		--from=markdown \
		--to=latex \
		--pdf-engine=xelatex \
		--include-in-header=template/si_header.tex \
		--bibliography=bib/references.bib \
		--csl=template/new-phytologist.csl \
		-o pdf/manuscript.pdf \
		ms/front_matter.md \
		ms/abstract.md \
		ms/introduction.md \
		ms/methods.md \
		ms/results.md \
		ms/discussion.md \
		ms/end_matter.md \
		template/ref_loc.md \
		ms/figure_legends.md \
		ms/si_list.md


pdf/supporting_information.pdf: ms/supplementary_figure_legends.md ms/supplementary_table_captions.md
	pandoc \
		--from=markdown \
		--to=latex \
		--pdf-engine=xelatex \
		--include-in-header=template/si_header.tex \
		--bibliography=bib/references.bib \
		--csl=template/new-phytologist.csl \
		-o pdf/supporting_information.pdf \
		ms/supplementary_figure_legends.md \
		ms/supplementary_table_captions.md

docx/supporting_information.docx: ms/supplementary_figure_legends.md ms/supplementary_table_captions.md
	pandoc \
		--from=markdown \
		--to=docx \
		--reference-doc=template/reference.docx \
		--bibliography=bib/references.bib \
		--csl=template/new-phytologist.csl \
		-o docx/supporting_information.docx \
		ms/supplementary_figure_legends.md \
		ms/supplementary_table_captions.md

