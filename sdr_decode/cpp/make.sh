#
#

target='loraphy'
c_files='main.cpp LoraPhy.cpp'
gnu_cc_ver='gnu++17'  # gnu++11, gnu++14, gnu++17, gnu++20, gnu++2b

rm -f "$target"
g++ $c_files -o "$target" -std="$gnu_cc_ver" -v

