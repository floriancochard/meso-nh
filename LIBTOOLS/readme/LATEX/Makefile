tools.ps : tools.dvi
	dvips -o $@ $<

tools.dvi : tools.tex
	latex tools.tex
	latex tools.tex
	latex tools.tex

clean:
	rm -f *.aux *.log *.toc *.dvi

realclean: clean
	rm -f *.ps
