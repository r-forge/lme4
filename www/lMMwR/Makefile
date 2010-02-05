R_PROG = R-devel
#R_PROG = R
FRONT = DMBmacros.tex acronym.tex dedication.tex preface.tex # acknow.tex foreword.tex
CHAPTERS = ChIntro ChMultiple ChCovariates ChLongitudinal ChTheory
APPENDICES =
RNW = $(CHAPTERS) $(APPENDICES)

## Currently unused, using R code to produce (and test for existance):
## Need some graphics as 'png' rather than 'pdf' to make the resulting pdf
## load *fast*
lrgGRAPHS = #figs/Intro-fm1profpair.pdf \
	    #figs/Multiple-fm2prpairs.pdf
smlGRAPHS= $(lrgGRAPHS:.pdf=.png)



all:  $(CHAPTERS:=.tex) $(RNW:=.R) $(RNW:=.tex) lMMwR.pdf lrgprt.pdf
# keep *.aux at least (for subsets in \includeonly{}):
#	${RM} *.out Rplots.ps figs/*.eps \
#		*.aux *.log *.toc *.nav *.snm *.vrb *.dvi

%.tex: %.Rnw
	${R_PROG} CMD Sweave $<

%.R: %.Rnw
	${R_PROG} CMD Stangle $<

%.pdf : %.tex
	${R_PROG} CMD texi2pdf $<

%.png: %.pdf
	convert $< $@


lMMwR.pdf: lMMwR.tex lMMwR.bib lMMwR.ind $(FRONT) $(RNW:=.tex) $(smlGRAPHS)

lrgprt.tex: lMMwR.pdf
	@echo "\\documentclass[letterpaper]{article}" > $@
	@echo "\\usepackage{pdfpages}" >> $@
	@echo "\\begin{document}" >> $@
	@echo "\\includepdf[nup=1x1,pages=-,noautoscale=true,scale=1.25,offset=0mm -12mm]{lMMwR}" >> $@
	@echo "\\end{document}" >> $@


lMMwR.ind: lMMwR.idx
	makeindex -c -s svind.ist lMMwR.idx

lMMwR.idx:
	touch lMMwR.idx


clean:
	$(RM) -f *.aux *.bbl *.blg *.idx *.ilg *.log *.ind *.toc *.lof *.lot
	$(RM) -f Rplots.* figs/*.pdf

realclean: clean
	$(RM) -f *~ $(RNW:=.tex) *.pdf *.rda lrgprt.tex figs/*.png

## TODO?  Use 'dist: clean' or 'dist: realclean'  and maybe *some* extra?
dist:
	$(RM) Rplots.ps $(RNW:=.R) $(CHAPTERS:=.tex) $(CHAPTERS:=.out) allR \
	   figs/*.pdf figs/*.eps *.aux *.log *.toc *.nav *.snm *.vrb *.dvi
