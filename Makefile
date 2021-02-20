SI.pdf: SI.tex refs.bib shareddefs.tex
	pdflatex $<
	bibtex $(subst .tex,,$<)
	pdflatex $<
	pdflatex $<

MT.pdf: MT.tex refs.bib shareddefs.tex auths.tex
	pdflatex $<
	bibtex $(subst .tex,,$<)
	pdflatex $<
	pdflatex $<

MT.docx: MT.md
	pandoc -s -o $@ $<

MT.md: MT.tex refs.bib shareddefs.tex auths.tex
	pandoc -s -o $@ --metadata bibliography=$(filter %.bib,$^) $<
