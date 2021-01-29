#!/usr/bin/env bash                                                                                                                                            
 
module load ruby
 
cmd="module list -l 2>&1 | grep -v default | tail -n +3 | cut -d \'/\' -f 1"

for mod in cdo python ncl; do
  #TODO echo ruby -e "mod='${mod}'; list = IO.popen('${cmd}').read.split; puts list.join('|')"
  foundMod=$(ruby -e "mod='${mod}'; list = IO.popen('${cmd}').read.split; puts list.include?(mod)")
  echo "Found $mod?: $foundMod"
  if [[ "true" = "${foundMod}" ]] ; then 
    echo "found '$mod' module"; 
  else 
    echo "could not find $mod"; 
    exit 1; 
  fi
done
exit 0
