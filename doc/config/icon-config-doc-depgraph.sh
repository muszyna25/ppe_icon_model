#!/bin/bash

script_dir=$(cd "$(dirname "$0")"; pwd)

in_file="${script_dir}/icon-config-doc-depgraph.dot"

out_format='svg'
out_file="${script_dir}/icon-config-doc-depgraph.${out_format}"

dot -Gconcentrate=true -T"${out_format}" "${in_file}" -o "${out_file}"
