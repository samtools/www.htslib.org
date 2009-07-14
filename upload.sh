./cleanup.sh
find . -name "*.html" | xargs ./csi2ssi.pl
mv index.html index.true.html
rsync -Cavz * lh3lh3,samtools@web.sourceforge.net:htdocs/
./cleanup.sh