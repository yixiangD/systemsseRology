format:
	R --slave -e 'styler::style_dir("R")'
	R --slave -e 'roxygen2::roxygenise()'
