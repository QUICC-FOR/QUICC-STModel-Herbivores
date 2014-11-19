currentState: draft.pdf

draft.pdf: tex/manuscript.tex tex/herbiveg.bib tex/herbiveg.bst tex/draft.tex graphs/*.pdf
	cd tex; pdflatex manuscript; bibtex manuscript; pdflatex manuscript; pdflatex manuscript 
	mv tex/manuscript.pdf draft.pdf

clean:
	cd tex; rm -f *.aux *.log *.bbl *.blg *.loc *.gz *.cut