#!/usr/bin/env ruby
require 'pp'
require "numru/netcdf"
require 'numru/netcdf_miss'
include NumRu

# 
# This script splits the glocal cell grid into northern and southern hemisphere.
#
# Requirements: ruby, ruby modules: ruby-netcdf

iFilename = ARGV[0]

iceVar    = 'p_ice_hi'
lon,lat   = 'clon','clat'

iFile     = NetCDF.open(iFilename)
lats      = iFile.var(lat).get
iceValues = iFile.var(iceVar).get_with_miss

# compute location indices and corresponding values
nhIndeces = (lats>0.0).where
shIndeces = (lats<0.0).where
iceValuesNH = iceValues[nhIndeces,true,1..-1]
iceValuesSH = iceValues[shIndeces,true,1..-1]


# create output
nhFile,shFile = "nh_#{iFilename}","sh_#{iFilename}"
[nhFile,shFile].each_with_index {|file,i|
  indeces = [nhIndeces,shIndeces][i]
  f = NetCDF.create(file)
  iFile.each_dim {|dim| 
    next if ['clon','clat','ncells'].include?(dim.name) or 
    f.def_dim(dim.name,dim.length)
  }
  ["clon","clat","ncells"].each {|hdim| 
    f.def_dim(hdim,indeces.size)
  }
 
  iFile.each_var{|var|
    a = {var.name => var.dim_names}
    pp a
    newvar = f.def_var( var.name, var.ntype, var.dim_names )
    var.each_att{|att| newvar.put_att( att.name, att.get )}
  }
  f.enddef
  iFile.each_var{|var| 
    puts var.name
    case var.name 
    when 'p_ice_hi'
    f.var(var.name).put(var.get[indeces,true,true])
    when 'p_ice_concSum'
    f.var(var.name).put(var.get[indeces,true])
    when 'cell_area'
    f.var(var.name).put(var.get[indeces])
    when 'clon','clat'
    f.var(var.name).put(var.get[indeces])
    when 'clon_vertices','clat_vertices'
    f.var(var.name).put(var.get[true,indeces])
    else
    f.var(var.name).put(var.get)
    end
  }
  f.close
}
