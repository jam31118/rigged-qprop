#!/bin/bash

### Parse command-line argument and set some configurations
# The first argument '$1' is absolute directory path from which archive will be copied
if [ -z "$1" ]
then echo "[ERROR] Please enter absolute directory path from which archive will be copied"; exit 1; 
else dirFrom="$1"
fi

# The second argument '$2' is absolute directory path to which archive will be copied
if [ -z "$2" ]
then echo "[ERROR] Please enter absolute directory path to which archive will be copied"; exit 1; 
else dirTo="$2"
fi

# The third argument '$3' is a whitespace-separated list of filename expression of files to exclude when copying archive
# .. The exclusion of file is required when the file is too large or needless etc.
excludeList="*.raw hydrogen_re-wf.dat"
if [ -n "$3" ]; then excludeList="$3"; fi


### Construct command
command="rsync -av $dirFrom $dirTo"
for x in $excludeList
do
	command=$command" "
	command=$command"--exclude=\"$x\""
done


### Run command
echo $command
eval $command
