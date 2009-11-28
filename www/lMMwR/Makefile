R_PROG = R-devel
#R_PROG = R
FRONT = DMBmacros.tex dedication.tex preface.tex
CHAPTERS = ChIntro ChMultiple ChLongitudinal
APPENDICES = 

all:  $(CHAPTERS:=.tex) $(CHAPTERS:=.R) $(APPENDICES:=.tex) $(APPENDICES:=.R) lMMwR.pdf lrgprt.pdf
	${RM} *.out Rplots.ps figs/*.eps \
          *.aux *.log *.toc *.nav *.snm *.vrb *.dvi

%.tex: %.Rnw
	${R_PROG} CMD Sweave $<

%.R: %.Rnw
	${R_PROG} CMD Stangle $<

%.pdf : %.tex
	${R_PROG} CMD texi2pdf $<

lMMwR.pdf: lMMwR.tex lMMwR.bib lMMwR.ind ${FRONT} ${CHAPTERS:=.tex} ${APPENDICES:=.tex}

lrgprt.tex: lMMwR.pdf
	@echo "\\documentclass[letterpaper]{article}" > $@
	@echo "\\usepackage{pdfpages}" >> $@
	@echo "\\\\begin{document}" >> $@
	@echo "\\includepdf[nup=1x1,pages=-,noautoscale=true,scale=1.25,offset=0mm -12mm]{lMMwR}" >> $@
	@echo "\\end{document}" >> $@


lMMwR.ind: lMMwR.idx
	makeindex -c -s svind.ist lMMwR.idx

lMMwR.idx:
	touch lMMwR.idx


clean:
	$(RM) -f *.aux *.bbl *.blg *.idx *.ilg *.pdf *.log *.ind *.toc *.lof *.lot
	$(RM) -f Rplots.ps *~ figs/*.pdf
	$(RM) $(CHAPTERS:=.tex)

dist:
	$(RM) Rplots.ps $(CHAPTERS:=.tex) $(CHAPTERS:=.out) allR \
	   figs/*.pdf figs/*.eps *.aux *.log *.toc *.nav *.snm *.vrb *.dvi
