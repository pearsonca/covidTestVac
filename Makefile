
R = Rscript $^ $@

grid.rds: gen_grid.R
	${R}

costs.png: fig_costs.R grid.rds
	${R}

figs: costs.png people.png protection.png

SI.pdf: SI.tex refs.bib shareddefs.tex auths.tex figs
	pdflatex $<
	bibtex $(subst .tex,,$<)
	pdflatex $<
	pdflatex $<