#1 - folder name

# should be input folder inside of folder containing ctg inputs.
mkdir ./output
for filename in $1/*.in; do
	echo "Running $filename"
	./count_tmp_graphlets "$filename" 6 4 1 "./output/$(basename "$filename" .in).out" -v ./vector_files/graphlets_6_4.txt -g ./vector_files/orbits_6_4.txt
	echo "Done with $filename"
done

