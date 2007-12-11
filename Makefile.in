EXECS=nph dispatcher plotmc1 plot_gixos gixos_nph
WEBROOT=/var/www/localhost

.PHONY: all

all:
	$(MAKE) -C cgi-bin

install:
	for f0 in $(EXECS); do \
		cp -v cgi-bin/$(f0) $(WEBROOT)/cgi-bin/; \
		done; \
	cp -v htdocs/*.php $(WEBROOT)/htdocs/; \
	cp -v htdocs/*.html $(WEBROOT)/htdocs/; \
	mkdir -p $(WEBROOT)/downloads

