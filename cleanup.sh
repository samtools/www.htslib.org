find . -name "*.shtml" | xargs rm
find . -name "*~" | xargs rm
if [ -f index.true.html ]; then
	mv -f index.true.html index.html
fi