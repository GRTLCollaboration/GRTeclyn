file(REMOVE_RECURSE
  "lib/libzmq.pdb"
  "lib/libzmq.so"
  "lib/libzmq.so.5"
  "lib/libzmq.so.5.2.5"
)

# Per-language clean rules from dependency scanning.
foreach(lang C CXX)
  include(CMakeFiles/libzmq.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
