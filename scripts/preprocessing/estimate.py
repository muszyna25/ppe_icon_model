#! /usr/bin/env python
#------------------------------------------------------------------------------
import getopt, sys, math

#"R<grid_refine_root>B<grid_refine_level>

#------------------------------------------------------------------------------
# set grid parameters
grid_refine_root  = 2
grid_refine_level = 3
vertical_layers   = 95

#------------------------------------------------------------------------------
# set machine parameters
#nodes_list=[ 32.0 ]
#opnemp_threads_list=[1.0, 4.0, 8.0, 16.0, 32.0]
#opnemp_threads_list=[ 8.0 ]
nproma_list=[ 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0,
34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0  ]

nodes_list=[   2.0 ]
openmp_threads_list=[ 4.0 ]
# 64=compute decompositions SMT mode
# 32=compute decompositions ST mode
total_threads_per_node=64.0

#------------------------------------------------------------------------------
# compute size of grid
cells = (20*grid_refine_root*grid_refine_root)*(4**grid_refine_level)
verts = (10*grid_refine_root*grid_refine_root)*4**grid_refine_level+2
edges = (30*grid_refine_root*grid_refine_root)*4**grid_refine_level

print "==============================================================="
print "= grid info: R",grid_refine_root,"B",grid_refine_level
print "  cells=",cells, " edges=", edges, " vertices=", verts
print "  vertical_layers=", vertical_layers
print "==============================================================="
#------------------------------------------------------------------------------

def display_thread_load(message, entities, threads, block_size):
  blocks = entities / block_size
  blocks_per_thread = blocks / threads
  blocks_per_thread_i = int(blocks_per_thread + 0.49)
  remain_block = (int(blocks + 0.49)) % (int(threads))
  #------------------------------------------------------------------------------
  print message, " blocks=",blocks, " blocks_per_thread=", blocks_per_thread,\
  blocks_per_thread_i, " remain_block=", remain_block

def compute_decomposition(nodes, opnemp_threads):
  mpi_procs_per_node = total_threads_per_node / opnemp_threads
  subdomains = mpi_procs_per_node * nodes
  owned_cells = cells / subdomains
  halo_cells = 4.4 * math.sqrt(owned_cells)
  domain_cells = owned_cells + halo_cells
  owned_edges = edges / subdomains
  owned_verts = verts / subdomains


  #------------------------------------------------------------------------------
  print "==============================================================="
  print "= decomposition info"
  print "  nodes=",nodes, " openmp threads=", opnemp_threads, " mpi_procs_per_node=", mpi_procs_per_node
  print "  subdomains=", subdomains
  print "  Per domain:  own cells=",owned_cells, " halo cells=", halo_cells, "total_cells=", domain_cells
  print "  Per domain:  own edges=", owned_edges, " own vertices=", owned_verts
  print "=========================="
  #------------------------------------------------------------------------------
  # compute nproma parameters
  for nproma in nproma_list:
    blocks_own_cell   = owned_cells   / nproma
    blocks_total_cell = domain_cells  / nproma
    blocks_own_edges  = owned_edges   / nproma
    blocks_own_verts  = owned_verts   / nproma
    #------------------------------------------------------------------------------
    print "= threads info, nproma=", nproma
    display_thread_load("own cells",   owned_cells, opnemp_threads, nproma)
    display_thread_load("total cells", domain_cells, opnemp_threads, nproma)
    display_thread_load("own edges",   owned_edges, opnemp_threads, nproma)
    display_thread_load("own verts",   owned_verts, opnemp_threads, nproma)
    print "=========================="

#------------------------------------------------------------------------------

for no_nodes in nodes_list:
  for omp_threads in openmp_threads_list:
    compute_decomposition(no_nodes, omp_threads)


sys.exit(0)


levels = 95
fields = 1

bytesize = 8

kilo = 10e6
mega = 10e9
giga = 10e12
tera = 10e15
peta = 10e18


size_per_3dvar = (gridpoints*levels*bytesize)/giga
size_per_step  = (gridpoints*fields*levels*bytesize)/giga
size_per_day = size_per_step * 4 # assuming output each minute
size_per_year = size_per_day * 365 # assuming output each minute

print "Model output estimate:"
print " - grid points per variable and level:     %10d" % gridpoints
print " - size per variable and horizontal slice: %10.3f gb" % size_per_3dvar
print " - size per output time:                   %10.3f gb" % size_per_step  
print " - size per output day:                    %10.3f gb" % size_per_day  
print " - size per output year:                   %10.3f gb" % size_per_year
