manuscript: docx/manuscript.docx 
si: pdf/supporting_information.pdf docx/supporting_information.docx
ms_pdf: pdf/manuscript.pdf
response: docx/response.docx

docx/manuscript.docx: ms/front_matter.md ms/abstract.md ms/introduction.md ms/methods.md ms/results.md ms/discussion.md ms/end_matter.md template/ref_loc.md ms/figure_legends.md ms/si_list.md bib/references.bib template/reference.docx template/journal-of-experimental-botany.csl template/ref_loc.md
	pandoc --from=markdown --to=markdown \
		ms/front_matter.md \
		ms/abstract.md \
		ms/introduction.md \
		ms/methods.md \
		ms/results.md \
		ms/discussion.md \
		ms/si_list.md \
		ms/end_matter.md \
		template/ref_loc.md \
		ms/figure_legends.md \
		| grep -vh "^!" \
		| pandoc --reference-doc=template/reference.docx \
		--filter pandoc-citeproc \
		--from=markdown \
		--to=docx \
		--bibliography=bib/references.bib \
		--csl=template/journal-of-experimental-botany.csl \
		-o docx/manuscript.docx

docx/response.docx: notes/reviewer_response.md bib/references.bib template/reference.docx template/journal-of-experimental-botany.csl template/ref_loc.md
	pandoc \
		--from=markdown \
		--to=docx \
		--reference-doc=template/reference.docx \
		-o docx/response.docx \
		notes/reviewer_response.md


pdf/manuscript.pdf: ms/front_matter.md ms/abstract.md ms/introduction.md ms/methods.md ms/results.md ms/discussion.md ms/end_matter.md template/ref_loc.md ms/supplementary_figure_legends.md bib/references.bib template/reference.docx template/journal-of-experimental-botany.csl ms/supplementary_table_captions.md template/ref_loc.md
	pandoc \
		--from=markdown \
		--to=latex \
		--pdf-engine=xelatex \
		--include-in-header=template/si_header.tex \
		--filter pandoc-include \
		--bibliography=bib/references.bib \
		--csl=template/journal-of-experimental-botany.csl \
		-o pdf/manuscript.pdf \
		ms/front_matter.md \
		ms/abstract.md \
		ms/introduction.md \
		ms/methods.md \
		ms/results.md \
		ms/discussion.md \
		ms/end_matter.md \
		template/ref_loc.md \
		ms/supplementary_figure_legends.md \
		ms/supplementary_table_captions.md


pdf/manuscript_mp.pdf: ms/abstract.md ms/introduction.md ms/methods.md ms/results.md ms/discussion.md template/ref_loc.md bib/references.bib template/reference.docx template/journal-of-experimental-botany.csl template/ref_loc.md ms/supplementary_figure_legends.md ms/supplementary_table_captions.md
	pandoc \
		--from=markdown \
		--to=latex \
		--pdf-engine=xelatex \
		--include-in-header=template/si_header.tex \
		--filter pandoc-include \
		--bibliography=bib/references.bib \
		--csl=template/journal-of-experimental-botany.csl \
		-o pdf/manuscript_mp.pdf \
		ms/abstract.md \
		ms/introduction.md \
		ms/methods.md \
		ms/results.md \
		ms/discussion.md \
		template/ref_loc.md \
		ms/supplementary_figure_legends.md \
		ms/supplementary_table_captions.md



pdf/supporting_information.pdf: ms/supplementary_figure_legends.md ms/supplementary_table_captions.md template/ref_loc.md
	pandoc \
		--from=markdown \
		--to=latex \
		--pdf-engine=xelatex \
		--include-in-header=template/si_header.tex \
		--bibliography=bib/references.bib \
		--csl=template/journal-of-experimental-botany.csl \
		-o pdf/supporting_information.pdf \
		ms/supplementary_figure_legends.md \
		ms/supplementary_table_captions.md \
		template/ref_loc.md


docx/supporting_information.docx: ms/supplementary_figure_legends.md ms/supplementary_table_captions.md template/ref_loc.md
	pandoc \
		--from=markdown \
		--to=docx \
		--reference-doc=template/reference.docx \
		--bibliography=bib/references.bib \
		--csl=template/journal-of-experimental-botany.csl \
		-o docx/supporting_information.docx \
		ms/supplementary_figure_legends.md \
		ms/supplementary_table_captions.md \
		template/ref_loc.md


newphyto_figs := pdf/newphyto/Figure_1_newphyto.pdf pdf/newphyto/Figure_2_newphyto.pdf pdf/newphyto/Figure_3_newphyto.pdf pdf/newphyto/Figure_4_newphyto.pdf pdf/newphyto/Figure_5_newphyto.pdf pdf/newphyto/Figure_6_newphyto.pdf
.PHONY: newphyto_figs
newphyto_figs: $(newphyto_figs)


 $(newphyto_figs): pdf/newphyto/Figure_%_newphyto.pdf: ms/newphyto_figs/figure_%_newphyto.md
	pandoc \
		--from=markdown \
		--to=latex \
		--pdf-engine=xelatex \
		--include-in-header=template/si_header.tex \
		--bibliography=bib/references.bib \
		--csl=template/new-phytologist.csl \
		-o $@ \
		$^
