#!/bin/bash
obj_prefix="tsurff-polar"
obj_suffix="dat"
obj="$obj_prefix.$obj_suffix"

if [ -f "$obj" ]; then echo "There is $obj alreadly."
echo "deleting $obj . . ."; rm "$obj"
else echo "Theres is no $obj."; fi

for f in "$obj_prefix"*
do
	echo "Combining \"$f\" to \"$obj\""
	cat "$f" >> "$obj"
done

