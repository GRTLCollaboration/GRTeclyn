#!/usr/bin/env bash

set -eu -o pipefail

find . -type d \( -name .git \
                  -o -name tmp_build_dir \
               \) -prune -o \
       -type f \( -name "*.f" -o -name "*.F" -o -name "*.f90" -o -name "*.F90" \
                  -o -name "*.py" \
                  -o -name "*.rst" \
                  -o -name "*.sh" \
                  -o -name "*.tex" \
                  -o -name "*.txt" \
                  -o -name "*.yml" \
               \) \
    -exec grep -Iq . {} \; \
    -exec sed -i 's/[[:blank:]]\+$//g' {} +

gitdiff=`git diff`

if [ -z "$gitdiff" ]
then
    exit 0
else
    echo -e "\nTrailing whitespaces at the end of a line are not allowed. Changes suggested by"
    echo -e "  ${0}\n"
    git --no-pager diff
    echo ""
    exit 1
fi
