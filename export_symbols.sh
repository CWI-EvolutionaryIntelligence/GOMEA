str=$(pip3 show gomea | grep "Location:")
libdir="${str/Location: /}/"
echo $libdir

#> symbols.txt
#> symbols-unmangled.txt
for f in $libdir/gomea/*.so
do
	echo $f #>> symbols.txt
	echo "Symbols"
	echo "============================================"
	nm $f #>> symbols.txt
	echo $f #>> symbols-unmangled.txt
	echo "Symbols unmangled"
	echo "============================================"
	nm -C $f #>> symbols-unmangled.txt
done
