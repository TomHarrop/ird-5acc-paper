################## SCRIPTS #####################

# Scripts that analyze raw-data
ANALYZE_SRC=$(wildcard analyze-data/*.R)
ANALYZE_SRC_IN_FOLDER=$(subst analyze-data/,,$(ANALYZE_SRC))

# Scripts that make figures
FIG_SRC=$(wildcard src/*.R)
FIG_SRC_IN_FOLDER=$(subst src/,,$(FIG_SRC))

############# DATA and FIGURES #################

# All tidy data
TIDY_DATA=$(wildcard data/*)

# All figures
FIGURES=$(wildcard fig/*)

# All raw daa
RAW_DATA=$(wildcard data-raw/*)

################## MANUSCRIPT  ##################

# dependencies
MAN_DEPENDS = \
  $(FIGURES) \
  template/reference.docx \
  bib/references.bib \
  template/genome-research.csl

# pandoc commands
DOCX_EXE = pandoc --reference-doc=template/reference.docx \
	--from=markdown \
	--to=docx \
	--bibliography=bib/references.bib \
	--output=$@ \


manuscript: docx/manuscript.docx 

docx/manuscript.docx: ms/front_matter.md ms/abstract.md ms/introduction.md ms/results.md ms/discussion.md ms/methods.md ms/figure_table_legends.md bib/references.bib template/reference.docx template/new-phytologist.csl template/ref_loc.md
	pandoc --reference-doc=template/reference.docx \
		--filter pandoc-citeproc \
		--from=markdown \
		--to=docx \
		--bibliography=bib/references.bib \
		--csl template/new-phytologist.csl \
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


################# MAKE FIGURES ####################

.PHONY: figures

figures:
	cd src; \
	for scr in $(FIG_SRC_IN_FOLDER); do \
	echo $$scr; \
	Rscript --vanilla $$scr; \
	done 

################ ANALYZE DATA ###################

.PHONY: analysis

analisys: 
	cd analyze-data; \
	for scr in $(ANALYZE_SRC_IN_FOLDER); do \
	echo $$scr; \
	Rscript --vanilla $$scr; \
	done 

