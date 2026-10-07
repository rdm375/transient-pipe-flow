.PHONY: docs clean

docs:
	cd docs && latexmk -pdf -interaction=nonstopmode -halt-on-error model.tex

clean:
	cd docs && latexmk -C model.tex
