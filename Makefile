DATA_RAW=$(wildcard data-raw/*)
ANALYZE_SRC=$(wildcard src/analyze*.R)
DATA_TIDY=$(patsubst src/analyze-%.R, data/%.Rdata, $(ANALYZE_SRC))
FIGS_SRC=$(wildcard src/fig*.R)
FIGS=$(patsubst src/%.R, fig/%.svg, $(FIGS_SRC))
TABS_SRC=$(wildcard src/table*.R)
TABS=$(patsubst src/%.R, tables/%.svg, $(TABS_SRC))

# For simplicity, specify common DATA dependecies for figures and tables
FIG_TAB_DEPEND = \
  $(DATA_TIDY) \
  src/helper-functions.R

# For simplicity, specify common dependencies for both the manuscript and the
# outline
TEXT_DEPEND = \
  $(FIGS) \
  $(TABS) \
  template/reference.docx \
  bib/references.bib \
  template/genome-research.csl \
  css/pandoc.css

# pandoc commands
DOCX_EXE = pandoc --reference-doc=template/reference.docx \
	--from=markdown \
	--to=docx \
	--bibliography=bib/references.bib \
	--output=$@ \

## all: make outline and manuscript
# the outline is just for us
all: outline manuscript

## outline: outline
outline: docx/outline.docx

docx/outline.docx: notes/outline.md $(TEXT_DEPEND)
	$(DOCX_EXE) $<

## manuscript:: make manuscript in docx and html
# The official manuscript is in docx, but the html version is much easier to use
.PHONY: manuscript
manuscript: docx/manuscript.docx html/manuscript.html

docx/manuscript.docx: ms/manuscript.md $(TEXT_DEPEND)
	$(DOCX_EXE) $<

html/manuscript.html: ms/manuscript.md $(TEXT_DEPEND)
	pandoc --from=markdown \
		--to=html \
		--toc\
		--css=../css/pandoc.css \
		--bibliography=bib/references.bib \
		--output=$@ \
		$<

## figs: read tidy data and make figures
.PHONY: figs
figs: $(FIGS)

fig/fig-%.svg: src/fig-%.R $(FIG_TAB_DEPEND)
	cd $(<D); Rscript --vanilla $(<F)

## tabs: read tidy data and make tables
.PHONY: tabs
tabs: $(TABS)

tables/table-%.svg: src/table-%.R $(FIG_TAB_DEPEND)
	cd $(<D); Rscript --vanilla $(<F)

## dats: analyze raw data and save tidy data inn Rdata format
#  Starting point of the analysis from which everything else depends
.PHONY: dats
dats: $(DATA_TIDY)

data/%.Rdata: src/analyze-%.R $(DATA-RAW)
	cd $(<D); Rscript --vanilla $(<F)

## variables   : Print variables.
.PHONY : variables
variables:
	@echo DATA-RAW: $(DATA_RAW)
	@echo ANALYZE_SRC: $(ANALYZE_SRC)
	@echo DATA_TIDY: $(DATA_TIDY)
	@echo FIGS_SRC: $(FIGS_SRC)
	@echo FIGS: $(FIGS)
	@echo TABS_SRC: $(TABS_SRC)
	@echo TABS: $(TABS)

## clean       : Remove auto-generated files.
.PHONY : clean
clean :
	rm -f $(DATA_TIDY)
	rm -f $(FIGS)
	rm -f $(TABS)
	rm -f html/manuscript.html

.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<
