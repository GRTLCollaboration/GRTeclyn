#!/bin/bash

#$1 = ${BINARY_DIR}
#$2 = ${CI_PROJECT_DIR}

count=0
while [ $count -le 9 ]; do
	if [[ -d plt00008 ]]; then
		"$1"/fcompare.gnu.ex --rel_tol 1e-10 --abs_tol 1e-10 --abort_if_not_all_found plt00008/ "${HOME}"/"$2"/.github/workflows/data/plt00008_compare
		"$1"/particle_compare.exe --rel_tol 1e-10 --abs_tol 1e-10 plt00008/ "${HOME}"/"$2"/.github/workflows/data/plt00008_compare/ particles
		break
	else
		((count++))
		echo "Attempt number " $count " of 10"
		sleep 30
	fi
done
