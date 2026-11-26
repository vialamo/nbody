# Makefile for generating docs

# Variables
CSS = source/style.css

BOOK_SRC = source/notes.md
BOOK_PDF = build/notes.pdf
BOOK_EPUB = build/notes.epub
BOOK_TITLE = "Notes on Cosmological Simulations"

METHODS_SRC = source/methods.md
METHODS_PDF = build/methods.pdf
METHODS_EPUB = build/methods.epub
METHODS_TITLE = "A Hybrid PM+PP/MFM Cosmological Hydrodynamics Code with GPU Acceleration"

.PHONY: all epub pdf book methods clean

# Default target
all: pdf

# Explicit targets
epub: $(BOOK_EPUB) $(METHODS_EPUB)
pdf: $(BOOK_PDF) $(METHODS_PDF)
book: $(BOOK_EPUB) $(BOOK_PDF)
methods: $(METHODS_EPUB) $(METHODS_PDF)

# Book targets
$(BOOK_EPUB): $(BOOK_SRC) $(CSS)
	@mkdir -p build
	pandoc $(BOOK_SRC) --css=$(CSS) -o $@ --mathml --to=epub3 -V documentclass=book --metadata title=$(BOOK_TITLE)

$(BOOK_PDF): $(BOOK_SRC)
	@mkdir -p build
	pandoc $(BOOK_SRC) -V geometry:"margin=0.8in" -V documentclass=book -V papersize=a4 --toc --toc-depth=3 -o $@

# Methods targets
$(METHODS_EPUB): $(METHODS_SRC) $(CSS)
	@mkdir -p build
	pandoc $(METHODS_SRC) -o $@ --mathml --to=epub3 -V documentclass=article --metadata title=$(METHODS_TITLE)

$(METHODS_PDF): $(METHODS_SRC)
	@mkdir -p build
	pandoc $(METHODS_SRC) -V geometry:"margin=0.8in" -V documentclass=article -V papersize=a4 -o $@

# Clean up generated files
clean:
	rm -rf build/
