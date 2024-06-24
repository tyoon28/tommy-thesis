#1 - folder name

# should be input folder inside of folder containing orca inputs.
mkdir ./output
for filename in $1/*.in; do
	echo "Running $filename"
	./orca 5 "$filename" "./output/$(basename "$filename" .in).out"
	echo "Done with $filename"
done

