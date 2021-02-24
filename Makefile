
R = Rscript $^ $@

main.png: fig_main.R
	${R}

figs: main.png

SI.pdf: SI.tex refs.bib shareddefs.tex auths.tex figs
	pdflatex $<
	bibtex $(subst .tex,,$<)
	pdflatex $<
	pdflatex $<