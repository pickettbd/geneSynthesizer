# MAKEFILE

# ------ VARIABLES -------------------
PYTHON3_BINARY=python
PREFIX=/usr/local/bin

# ------ TARGETS ---------------------
.PHONY: all clean realclean prep permissions install uninstall src/geneSynthesizer.py

all: prep bin/geneSynthesizer permissions

clean:
	@rm -rf bin/* || true
	@if [ -e $(PREFIX)/geneSynthesizer -a -r $(PREFIX)/geneSynthesizer -a -w $(PREFIX)/geneSynthesizer ] ; then rm -f $(PREFIX)/geneSynthesizer || true; fi
	@chmod 640 src/* || true
	@chmod 770 src/__pycache__ || true

realclean: clean
	@rm -rf bin src/__pycache__ src/*.pyc || true
	@chmod 640 src/*.py || true

prep:
	@mkdir -p bin || true

src/geneSynthesizer.py:
	@$(PYTHON3_BINARY) -c 'import py_compile as pc; pc.compile("'$@'")'

bin/geneSynthesizer: src/geneSynthesizer.py
	@ln -s ../$< $@ 2> /dev/null || true

permissions:
	@chmod 640 src/* || true
	@chmod 770 src/__pycache__ &> /dev/null || true
	@chmod 750 src/geneSynthesizer.py || true

install:
	@if [ -e bin/geneSynthesizer ]; then cp bin/geneSynthesizer $(PREFIX) || true; else echo "ERROR: \`bin/geneSynthesizer*' does not exist. Did you forget to run \`make' first?"; fi

uninstall:
	@if [ -e $(PREFIX)/geneSynthesizer ]; then rm $(PREFIX)/geneSynthesizer; fi
