.PHONY: clean doc.pdf html

doc.pdf:
	latexmk -pdf doc.tex

html:
	jupyter-book build --html

clean:
	rm -f doc.aux doc.bbl doc.blg doc.fdb_latexmk doc.pdf doc.synctex.gz doc.fls doc.idx doc.log doc.out doc.toc doc.xdv doc.ind doc.ilg
