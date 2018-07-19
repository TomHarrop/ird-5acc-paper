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


manuscrpit: docx/manuscript.docx 

docx/manuscript.docx: ms/manuscript.md $(MAN_DEPEND)
	$(DOCX_EXE) $<


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

