.PHONY: doc.pdf doc.html

doc.pdf:
	latexmk -pdf doc.tex

doc.html:
	pandoc doc.tex -o doc.html \
	  --metadata title="ALF Documentaion" \
	  --mathjax \
	  -H header-mathjax.html \
	  --standalone \
	  --toc \
	  --css=pandoc.css \
	  --number-sections
