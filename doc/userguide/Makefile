SHELL = /bin/ksh

USERG_BIBLIOGRAPHY = icon_user_guide_bib.bib
USERG_SOURCES      = icon_user_guide.tex authorlist.tex content/*.tex \
	             ../icon_nml_tables.tex ../icon_ocean_nml_tables.tex \
                     ./*.bib
USERG_STYLE        = style/chicago_lmdocu.bst style/icon_docu_style.cls \
	             style/pageformatting_toc.tex
USERG_PICTURES     = pictures/*png

all: icon_user_guide

icon_user_guide: icon_user_guide.pdf

icon_user_guide.pdf: $(USERG_SOURCES) $(USERG_STYLE) \
                     $(USERG_PICTURES) $(USERG_BIBLIOGRAPHY)
	pdflatex icon_user_guide
	bibtex icon_user_guide
	makeindex icon_user_guide
	pdflatex icon_user_guide
	pdflatex icon_user_guide
