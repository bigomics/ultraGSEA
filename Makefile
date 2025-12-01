build: doc vignettes
	R -e "devtools::build()"

doc:R/*.R
	R -e "devtools::document()"

vignettes: vignettes/*.Rmd
	R -e "devtools::build_vignettes()"

install: 
	R CMD INSTALL .

check: clean
	R -e 'devtools::check()'

biocheck: clean
	R -e 'BiocCheck::BiocCheck(".")'

test:
	R -e "devtools::test()"

clean:
	rm -f `find . -name '.\#*' -o -name '\#*' -o -name '*~' -printf '"%p" '`

echo:
	echo `find . -name '.\#*' -o -name '\#*' -o -name '*~' -printf '"%p" '`

FORCE: ;
