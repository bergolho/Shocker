# --------------------------------------------------------------------------------
# Memory leak test using Valgrind
# --------------------------------------------------------------------------------

#!/bin/bash
PNAME="./bin/Purkinje-Parser"
INPUT_FILE="inputs/simple_purkinje_network_1.vtk"
OUTPUT_FILE="outputs/simple_purkinje_network_1.pkje"

if [ ! -f $PNAME ]; then
	./recompile_project.sh
fi

valgrind --leak-check=full --show-leak-kinds=all ./$PNAME $INPUT_FILE $OUTPUT_FILE