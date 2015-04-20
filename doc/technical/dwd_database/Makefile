all:
	pdflatex database_main.tex
	bibtex database_main
	pdflatex database_main.tex
	pdflatex database_main.tex

clean:
	rm -f *.aux *.log *.out *.toc *~ *.bbl *.blg *.ver *.hst
