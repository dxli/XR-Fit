htdir = /var/www/localhost/htdocs
all:
	echo "Nothing to do"

install:
	for f0 in xray.html config.php post.php last.php gixos.html \
		post-gixos.php gixos-last.php ;\
	do install -m644 -o apache -d ${htdir} ${f0}; \
	done


