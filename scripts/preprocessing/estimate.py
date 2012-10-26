#! /usr/bin/env python
#------------------------------------------------------------------------------
import getopt, sys, math

#"R<grid_refine_root>B<grid_refine_level>

#------------------------------------------------------------------------------

cells=0
verts=0
edges=0

#------------------------------------------------------------------------------
vertical_levels=96
#------------------------------------------------------------------------------
# set machine parameters
#nodes_list=[ 32.0 ]
#opnemp_threads_list=[1.0, 4.0, 8.0, 16.0, 32.0]
#opnemp_threads_list=[ 8.0 ]
nproma_list=[  8.0, 10.0, 12.0, 16.0, 18.0 ]

#mpi_procs_list=[ 128.0, 256.0, 512.0, 541.0, 640.0, 768.0,  1082.0 ]
#mpi_nodes_list=[ 2.0, 4.0, 8.0, 16.0, 24.0, 34.0 ]
omp_threads=2.0
mpi_procs_pernode=32
# 64=compute decompositions SMT mode
# 32=compute decompositions ST mode
total_threads_per_node=64.0


def display_thread_load(message, entities, threads, block_size):
    #blocks_own_cell   = round(owned_cells/nproma)
    #blocks_own_cell   = (owned_cells/nproma)
    #blocks_total_cell = round(domain_cells  / nproma)
    #blocks_own_edges  = round(owned_edges   / nproma)
    #blocks_own_verts  = round(owned_verts   / nproma)
  
  blocks = (entities / block_size )
  blocks_i = round(blocks + 0.25)
  blocks_per_thread = blocks_i / threads
  blocks_per_thread_i = round(blocks_per_thread)
  iterations = int(entities / blocks_per_thread_i)
  #remain_block = entities - (blocks_per_thread_i * threads * iterations)
  remain_block = blocks_i % threads
  #------------------------------------------------------------------------------
  print message, " blocks=",blocks, " blocks_per_thread=", blocks_per_thread,\
  blocks_per_thread_i, " remain_block=", remain_block

def compute_decomposition(mpi_procs, opnemp_threads):
  subdomains = mpi_procs
  owned_cells = round((cells / subdomains) + 0.499999)
  halo_cells = 4.4 * math.sqrt(owned_cells)
  domain_cells = owned_cells + halo_cells
  owned_edges = edges / subdomains
  owned_verts = verts / subdomains


  #------------------------------------------------------------------------------
  print "==============================================================="
  print "= decomposition info"
  print "  mpi_procs=",mpi_procs, " openmp threads=", opnemp_threads
  print "  subdomains=", subdomains
  print "  Per domain:  own cells=",owned_cells, " halo cells=", halo_cells, "total_cells=", domain_cells
  print "  Per domain:  own edges=", owned_edges, " own vertices=", owned_verts
  print "=========================="
  #------------------------------------------------------------------------------
  # compute nproma parameters
  for nproma in nproma_list:
    #------------------------------------------------------------------------------
    print "= threads info, nproma=", nproma
    display_thread_load("own cells",   owned_cells, opnemp_threads, nproma)
    display_thread_load("total cells", domain_cells, opnemp_threads, nproma)
    display_thread_load("own edges",   owned_edges, opnemp_threads, nproma)
    display_thread_load("own verts",   owned_verts, opnemp_threads, nproma)
    print "=========================="

#------------------------------------------------------------------------------

def check_nproma(grid_name):

  print "==============================================================="
  print "= "
  print "= grid_name:",grid_name
  print "  cells=",cells, " edges=", edges, " vertices=", verts
  print "  vertical_levels=", vertical_levels
  print "= "
  print "==============================================================="
  #------------------------------------------------------------------------------
  #for mpi_procs in mpi_procs_list:
      #compute_decomposition(mpi_procs, omp_threads)
  for mpi_nodes in mpi_nodes_list:
      compute_decomposition(mpi_nodes * mpi_procs_pernode, omp_threads)


#------------------------------------------------------------------------------
# compute size of grid
def compute_icon_standard_size(grid_refine_root, grid_refine_level):
  cells = (20*grid_refine_root*grid_refine_root)*(4**grid_refine_level)
  verts = (10*grid_refine_root*grid_refine_root)*4**grid_refine_level+2
  edges = (30*grid_refine_root*grid_refine_root)*4**grid_refine_level
  return (cells, verts, edges)

mpi_nodes_list=[ 214.0 ]
(cells, verts, edges)=compute_icon_standard_size(2, 6)
check_nproma("iconR2B04-grid")

sys.exit(0)

mpi_nodes_list=[ 1.0, 4.0, 8.0, 16.0, 20.0, 27.0 ]
#nproma_list=[   16   12   20   10     8    12   ]
(cells, verts, edges)=compute_icon_standard_size(2, 4)
#check_nproma("iconR2B04-grid")

# iconR2B05-grid_dec-1082.nc
mpi_nodes_list=[ 1.0, 7.0, 12.0, 23.0 ]
#nproma_list=  [  12  12   12     9   ]
cells = 25920
verts = 12962
edges = 38880
#check_nproma("iconR2B05-grid_dec-1082")

# iconR2B05-grid_dec-362.nc
mpi_nodes_list=[ 1.0, 6.0, 27.0, 45.0 ]
#nproma_list=  [  12  18   10     12 ]
cells = 34560
verts = 17282
edges = 51840
#check_nproma("iconR2B05-grid_dec-362")

# iconR2B05-grid_dec-1922
mpi_nodes_list=[ 1.0, 6.0, 12.0,  24.0 ]
#nproma_list=  [  12  12    12    10   ]
cells = 46080
verts = 23042
edges = 69120
#check_nproma("iconR2B05-grid_dec-1922")
sys.exit(0)


#mpi_nodes_list=[  8.0 ]
## iconR2B06-grid_dec-1082.nc
#cells = 103680 
#verts = 51842
#edges = 155520
#check_nproma("iconR2B06-grid_dec-1082")
##nproma = [ 16 ]

## iconR2B06-grid_dec-1442.nc
#cells = 138240
#verts = 69122
#edges = 207360
#check_nproma("iconR2B06-grid_dec-1442")
##nproma = [ 12 ]

## iconR2B06-grid_dec-1922.nc
#cells = 184320
#verts = 92162
#edges = 276480
#check_nproma("iconR2B06-grid_dec-1922")
##nproma = [ 16 ]

## iconR2B07-grid_dec-2432.nc
#cells = 233280
#verts = 116642
#edges = 349920
#check_nproma("iconR2B07-grid_dec-2432")
##nproma = [ 16 ]


## iconR2B06-grid.nc
(cells, verts, edges)=compute_icon_standard_size(2, 6)
mpi_nodes_list=[  24.0 ]
check_nproma("iconR2B06-grid")
#nproma = [  16 ]

# iconR2B07_dec-4322.nc
#cells = 414720
#verts = 207362
#edges = 622080
#mpi_nodes_list=[  24.0 ]
#nproma = [  12 ]
#check_nproma("iconR2B07_dec-4322.nc")

sys.exit(0)

gridpoints=0.0
(cells, verts, edges)=compute_icon_standard_size(2, 6)
gridpoints=round(cells)

levels = 96.0
fields = 1.0

bytesize = 8.0

kilo = 1024.0
mega = kilo * kilo
giga = kilo * mega
tera = kilo * giga
peta = kilo * tera


size_per_3dvar = (gridpoints*levels*bytesize)/giga
size_per_step  = (gridpoints*fields*levels*bytesize)/giga
size_per_day = size_per_step * 4.0 # assuming output every six hours
size_per_year = size_per_day * 365.0 # 

print "Model output estimate:"
print " - grid points per variable and level:     %10d" % gridpoints
print " - size per variable and horizontal slice: %10.3f gb" % size_per_3dvar
print " - size per output time:                   %10.3f gb" % size_per_step  
print " - size per output day:                    %10.3f gb" % size_per_day  
print " - size per output year:                   %10.3f gb" % size_per_year
