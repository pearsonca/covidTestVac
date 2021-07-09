
main.pdf: WOR.tex refs.bib auths.tex main.png WellcomeOR_styles.sty WellcomeOR_logo_black.pdf
	pdflatex $<
	bibtex $(subst .tex,,$<)
	pdflatex $<
	pdflatex $<

R = Rscript $^ $@

main.png: fig_main.R
	${R}