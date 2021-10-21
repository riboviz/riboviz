## Makefile to create .svg files from .dot files using Graphviz 'dot'

DOT_DIR=docs/images
DOT_FILES=$(wildcard $(DOT_DIR)/*.dot)
SVG_DIR=docs/images
SVG_FILES=$(patsubst $(DOT_DIR)/%.dot, $(SVG_DIR)/%.svg, $(DOT_FILES))

$(SVG_DIR)/%.svg : $(DOT_DIR)/%.dot
	dot -Tsvg $< > $@

.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<

## clean       : Remove auto-generated files.
.PHONY : clean
clean :
	rm -f $(SVG_FILES)

## svgs        : Make SVG files from DOT files.
.PHONY : svgs
svgs : $(SVG_FILES)
