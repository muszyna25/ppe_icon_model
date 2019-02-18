/\s*SRCS\s*=/,/^\s*$/{//!d;s/^\(\s*SRCS\s*=\).*$/\1 \\/;/SRCS/a\
OCEFILES_DUMMY
}
/\s*OBJS\s*=/,/^\s*$/{//!d;s/^\(\s*OBJS\s*=\).*$/\1 $(patsubst %.f90,%.o,$(SRCS)) /; }
s,../bin/grid_command,,
s,../bin/test_divide_cell_mpi,,
s,../bin/test_divide_cell,,
/^\s*grid_command.o:/,/^\s*$/d
/^\s*test_divide_cell_mpi.o:/,/^\s*$/d
/^\s*test_divide_cell.o:/,/^\s*$/d
/^\s*clean\s*:.*/,${1p;d}

