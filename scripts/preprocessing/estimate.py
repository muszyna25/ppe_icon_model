#! /usr/bin/env python
#------------------------------------------------------------------------------
import getopt, sys, math

#"R<grid_refine_root>B<grid_refine_level>

#------------------------------------------------------------------------------
# set grid parameters
grid_refine_root  = 2
grid_refine_level = 4
vertical_layers   = 95

#------------------------------------------------------------------------------
# set machine parameters
#nodes_list=[ 32.0 ]
#opnemp_threads_list=[1.0, 4.0, 8.0, 16.0, 32.0]
#opnemp_threads_list=[ 8.0 ]
nproma_list=[ 6.0, 8.0, 12.0, 16.0, 18.0, 22.0, 24.0,  32.0 ]

mpi_procs_list=[ 128.0, 256.0, 512.0, 541.0, 640.0, 768.0,  1082.0 ]
mpi_nodes_list=[ 2.0, 4.0, 8.0, 16.0, 24.0, 34.0 ]
omp_threads=2.0
mpi_procs_pernode=32
# 64=compute decompositions SMT mode
# 32=compute decompositions ST mode
total_threads_per_node=64.0

#------------------------------------------------------------------------------
# compute size of grid
cells = (20*grid_refine_root*grid_refine_root)*(4**grid_refine_level)
verts = (10*grid_refine_root*grid_refine_root)*4**grid_refine_level+2
edges = (30*grid_refine_root*grid_refine_root)*4**grid_refine_level

# iconR2B05-grid_dec-1082.nc
#cells = 25920
#verts = 12962
#edges = 38880
# iconR2B05-grid_dec-362.nc
cells = 34560
verts = 17282
edges = 51840

print "==============================================================="
print "= grid info: R",grid_refine_root,"B",grid_refine_level
print "  cells=",cells, " edges=", edges, " vertices=", verts
print "  vertical_layers=", vertical_layers
print "==============================================================="
#------------------------------------------------------------------------------

def display_thread_load(message, entities, threads, block_size):
  blocks = round(entities / block_size+ 0.5)  
  blocks_per_thread = blocks / threads
  blocks_per_thread_i = round(blocks_per_thread)
  iterations = int(entities / blocks_per_thread_i)
  remain_block = entities - (blocks_per_thread_i * threads * iterations)
  remain_block = blocks % threads
  #------------------------------------------------------------------------------
  print message, " blocks=",blocks, " blocks_per_thread=", blocks_per_thread,\
  blocks_per_thread_i, " remain_block=", remain_block

def compute_decomposition(mpi_procs, opnemp_threads):
  subdomains = mpi_procs
  owned_cells = cells / subdomains
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
    blocks_own_cell   = round(owned_cells/nproma + 0.5)
    blocks_total_cell = round(domain_cells  / nproma+ 0.5)
    blocks_own_edges  = round(owned_edges   / nproma+ 0.5)
    blocks_own_verts  = round(owned_verts   / nproma+ 0.5)
    #------------------------------------------------------------------------------
    print "= threads info, nproma=", nproma
    display_thread_load("own cells",   owned_cells, opnemp_threads, nproma)
    display_thread_load("total cells", domain_cells, opnemp_threads, nproma)
    display_thread_load("own edges",   owned_edges, opnemp_threads, nproma)
    display_thread_load("own verts",   owned_verts, opnemp_threads, nproma)
    print "=========================="

#------------------------------------------------------------------------------

#for mpi_procs in mpi_procs_list:
    #compute_decomposition(mpi_procs, omp_threads)
for mpi_nodes in mpi_nodes_list:
    compute_decomposition(mpi_nodes * mpi_procs_pernode, omp_threads)


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
