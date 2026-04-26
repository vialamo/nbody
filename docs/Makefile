# Makefile for generating Cosmology Notes

# Variables
SRC = source/book.md
CSS = source/style.css
PDF_OUT = build/book.pdf
EPUB_OUT = build/book.epub
TITLE = "Notes on N-Body/Hydrodynamics Simulation"

# Default target runs both
all: epub pdf

# Explicit targets
epub: $(EPUB_OUT)
pdf: $(PDF_OUT)

# How to build the EPUB
$(EPUB_OUT): $(SRC) $(CSS)
	@mkdir -p build
	pandoc $(SRC) --css=$(CSS) -o $@ --mathml --to=epub3 -V documentclass=book --metadata title=$(TITLE)

# How to build the PDF
$(PDF_OUT): $(SRC)
	@mkdir -p build
	pandoc $(SRC) -V geometry:"margin=0.8in" -V documentclass=book -V papersize=a4 --toc --toc-depth=3 -o $@

# Clean up generated files
clean:
	rm -rf build/