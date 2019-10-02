#! /usr/bin/env perl
#__________________________________________________________________________________________________________________________________
#
# Creates Makefiles for the list of given source code directories. 
# This program is highly specialized for ICON and cannot be used 
# with other packages.
#__________________________________________________________________________________________________________________________________
#
use strict;
use warnings;
#
use Cwd;
use File::Copy;
use Getopt::Long;
use File::Path;
use File::Basename;
use File::Spec;
#__________________________________________________________________________________________________________________________________
# Option processing

my $target;
my $srcdirs;
my $enable_atmo;
my $enable_ocean;
my $enable_jsbach;
my $enable_rte_rrtmgp;
my $enable_testbed;
my $enable_serialization;
my $iconsrcdir = "src";
my $iconppdir;

our $enable_gpu;

GetOptions( 
	    'target=s'  => \$target,
            'srcdirs=s' => \$srcdirs,
            'enable_atmo=s' => \$enable_atmo,
            'enable_ocean=s' => \$enable_ocean,
            'enable_jsbach=s' => \$enable_jsbach,
            'enable_rte_rrtmgp=s' => \$enable_rte_rrtmgp,
            'enable_testbed=s' => \$enable_testbed,
            'enable_gpu=s' => \$enable_gpu, 
            'enable_serialization=s' => \$enable_serialization,
            'iconsrcdir=s' => \$iconsrcdir,
            'iconppdir=s' => \$iconppdir,
	    ) or die "\n\nUsage: config/createMakefiles.pl --target=<OS_CPU> --srcdirs=< list of src directories>\n";  

#__________________________________________________________________________________________________________________________________
#

$target =~ s/\s+//g;

$srcdirs =~ s/[\"\']*//g;
$srcdirs =~ s/\s+/ /g;
$srcdirs =~ s/^\s+//;

$iconsrcdir =~ s/\s+//g;

my (@directories) = split / /, $srcdirs;

#__________________________________________________________________________________________________________________________________
# determine base directory

my $prefix = &cwd;

#__________________________________________________________________________________________________________________________________
# make architecture dependend build directories

print "\n";
print "createMakefiles:\n\n";

my $build_path = &BuildSetup ($prefix, $target, \@directories);

if (defined $iconppdir) {
  $iconppdir = $prefix . "/" . $iconppdir;
} else {
  $iconppdir = $build_path . "/" . $iconsrcdir;
}

#__________________________________________________________________________________________________________________________________
# collect src files og rte-rrtmgp dependent on the accelerator type

if ( -d "externals/rte-rrtmgp") {

    # 0. select kernel type cpu(default), openacc
    my $kernel_type = $enable_rte_rrtmgp;

    if ($kernel_type ne "none") {
	print "\nProcess rte/rrtmgp for ${kernel_type} ...\n\n";

	# 1. create src directory
	if ( ! -d "externals/rte-rrtmgp/src" ) {
	    mkdir "externals/rte-rrtmgp/src";
	}
	# 2. copy all files from rte, rrtmgp, extensions to src
	opendir(DIR, "externals/rte-rrtmgp/rte");
	my @f90s = grep /\.(f90|F90)/, readdir(DIR);
	closedir(DIR);
	foreach my $f90 ( @f90s ) {
	    my $in = $f90;
	    my $out = $in;
	    $out =~ s/F90/f90/;
	    copy ( "externals/rte-rrtmgp/rte/${in}", "externals/rte-rrtmgp/src/${out}" );
	}
	opendir(DIR, "externals/rte-rrtmgp/rrtmgp");
	@f90s = grep /\.(f90|F90)/, readdir(DIR);
	closedir(DIR);
	foreach my $f90 ( @f90s ) {
	    my $in = $f90;
	    my $out = $in;
	    $out =~ s/F90/f90/;
	    copy ( "externals/rte-rrtmgp/rrtmgp/${in}", "externals/rte-rrtmgp/src/${out}" );
	}
	opendir(DIR, "externals/rte-rrtmgp/extensions");
	@f90s = grep /\.(f90|F90)/, readdir(DIR);
	closedir(DIR);
	foreach my $f90 ( @f90s ) {
	    my $in = $f90;
	    my $out = $in;
	    $out =~ s/F90/f90/;
	    copy ( "externals/rte-rrtmgp/extensions/${in}", "externals/rte-rrtmgp/src/${out}" );
	}
	opendir(DIR, "externals/rte-rrtmgp/extensions/cloud_optics");
	@f90s = grep /\.(f90|F90)/, readdir(DIR);
	closedir(DIR);
	foreach my $f90 ( @f90s ) {
	    my $in = $f90;
	    my $out = $in;
	    $out =~ s/F90/f90/;
	    copy ( "externals/rte-rrtmgp/extensions/cloud_optics/${in}", "externals/rte-rrtmgp/src/${out}" );
	}
	# 3. copy architecture related files
	if ($kernel_type eq "openacc") {
            opendir(DIR, "externals/rte-rrtmgp/rte/kernels-openacc");
            my @f90s = grep /\.(f90|F90)/, readdir(DIR);
            closedir(DIR);
            foreach my $f90 ( @f90s ) {
                my $in = $f90;
                my $out = $in;
                $out =~ s/F90/f90/;
                copy ( "externals/rte-rrtmgp/rte/kernels-openacc/${in}", "externals/rte-rrtmgp/src/${out}" );
            }
            opendir(DIR, "externals/rte-rrtmgp/rrtmgp/kernels-openacc");
            @f90s = grep /\.(f90|F90)/, readdir(DIR);
            closedir(DIR);
            foreach my $f90 ( @f90s ) {
                my $in = $f90;
                my $out = $in;
                $out =~ s/F90/f90/;
                copy ( "externals/rte-rrtmgp/rrtmgp/kernels-openacc/${in}", "externals/rte-rrtmgp/src/${out}" );
            }
        }
        else {
            opendir(DIR, "externals/rte-rrtmgp/rte/kernels");
            my @f90s = grep /\.(f90|F90)/, readdir(DIR);
            closedir(DIR);
            foreach my $f90 ( @f90s ) {
                my $in = $f90;
                my $out = $in;
                $out =~ s/F90/f90/;
                copy ( "externals/rte-rrtmgp/rte/kernels/${in}", "externals/rte-rrtmgp/src/${out}" );
            }
            opendir(DIR, "externals/rte-rrtmgp/rrtmgp/kernels");
            @f90s = grep /\.(f90|F90)/, readdir(DIR);
            closedir(DIR);
            foreach my $f90 ( @f90s ) {
                my $in = $f90;
                my $out = $in;
                $out =~ s/F90/f90/;
                copy ( "externals/rte-rrtmgp/rrtmgp/kernels/${in}", "externals/rte-rrtmgp/src/${out}" );
            }
        }
    }
}

#__________________________________________________________________________________________________________________________________
# collect include files in build tree

copy ("config/config.h", "${build_path}/include/config.h");

my @incs = "";

# as only cdi.inc is left and that has been replaced by the
# iso-c-binding module file it is not required anymore.
#
#opendir(DIR, "include") or die "Unable to open include:$!\n";
#my @incs = grep /\.(inc|h)/, readdir(DIR);
#closedir(DIR);
#foreach my $inc ( @incs ) {
#    copy ( "include/${inc}", "${build_path}/include/${inc}" );
#}


if ( -d "externals/mtime/include" ) {
    opendir(DIR, "externals/mtime/include");
    @incs = grep /\.(inc|h)/, readdir(DIR);
    closedir(DIR);
    foreach my $inc ( @incs ) {
	copy ( "externals/mtime/include/${inc}", "${build_path}/include/${inc}" );
    }
}

if ( -d "externals/ecrad/include" ) {
    opendir(DIR, "externals/ecrad/include");
    @incs = grep /\.(inc|h)/, readdir(DIR);
    closedir(DIR);
    foreach my $inc ( @incs ) {
        copy ( "externals/ecrad/include/${inc}", "${build_path}/include/${inc}" );
    }
}

if ( -d "externals/ecrad/radiation" ) {
    opendir(DIR, "externals/ecrad/radiation");
    @incs = grep /\.(inc|h)/, readdir(DIR);
    closedir(DIR);
    foreach my $inc ( @incs ) {
        copy ( "externals/ecrad/radiation/${inc}", "${build_path}/include/${inc}" );
    }
}

if ( -d "externals/tixi/include" ) {
    opendir(DIR, "externals/tixi/include");
    @incs = grep /\.(inc|h)/, readdir(DIR);
    closedir(DIR);
    foreach my $inc ( @incs ) {
	copy ( "externals/tixi/include/${inc}", "${build_path}/include/${inc}" );
    }
}

if ( -d "externals/yac/include" ) {
    opendir(DIR, "externals/yac/include");
    @incs = grep /\.(inc|h)/, readdir(DIR);
    closedir(DIR);
    foreach my $inc ( @incs ) {
	copy ( "externals/yac/include/${inc}", "${build_path}/include/${inc}" );
    }
}

if ( ($enable_jsbach eq "yes")
     and ! -d "$iconsrcdir/lnd_phy_jsbach"
     and -d 'externals/jsbach/src' ) {
    symlink File::Spec->abs2rel("$prefix/externals/jsbach/src", "$iconsrcdir"), "$iconsrcdir/lnd_phy_jsbach"; 
}

if ( ($enable_gpu eq "yes") and -d 'externals/cub' ) {
    symlink "$prefix/externals/cub", "${build_path}/include/cub";
}

if ( ($enable_ocean eq "yes") and -d "src/ocean/include" ) {
    opendir(DIR, "src/ocean/include");
    @incs = grep /\.(inc|h)/, readdir(DIR);
    closedir(DIR);
    foreach my $inc ( @incs ) {
	copy ( "src/ocean/include/${inc}", "${build_path}/include/${inc}" );
    }
}

#__________________________________________________________________________________________________________________________________
# scan dependencies (recursive)

my %vpath_directories  = ();
    
my %target_programs    = ();
    
my %module_definitions = ();
my %module_usage       = ();
    
my %fortran_includes   = ();
my %c_includes         = (); 
    
my %ifdefs             = ();
    
my @source_files       = ();

foreach my $dir ( @directories ) {

# global variables

    %vpath_directories  = ();
    
    %target_programs    = ();
    
    %module_definitions = ();
    %module_usage       = ();
    
    %fortran_includes   = ();
    %c_includes         = (); 
    
    %ifdefs             = ();
    
    @source_files       = ();

    &ScanDirectory ($dir, $dir, 0);

    my $print_path = $build_path;
    $print_path =~ s/$prefix\///;
    $dir =~ s/$print_path\///;

    my $iconrelsrcdir = $iconsrcdir;
    $iconrelsrcdir =~ s/$print_path\///;

    print "creating $print_path/$dir/Makefile\n";

    mkpath("$build_path/$dir");
    open MAKEFILE, ">$build_path/$dir/Makefile";

    print MAKEFILE <<EOF
#----------------------------------------------------------
# Makefile automatically generated by createMakefiles.pl
#----------------------------------------------------------

EOF
        ;

    my $add_vpath_level = 0;
    if ( "$dir" ne "$iconrelsrcdir" ) {
	if ( $dir =~ m/^externals/) {
	    my @subdirs = split(/\//, $dir);
	    print MAKEFILE "SHELL = /bin/bash\n\n",
                "LIB  = $subdirs[1]\n\n";
	    $add_vpath_level = 2;

	} else {
	    print MAKEFILE "LIB  = $dir\n\n";
	}
    }

#__________________________________________________________________________________________________________________________________
# write VPATH
    
    my @vpath = ();
    push @vpath, "VPATH = ";
    while ( my ($key, $value) = each(%vpath_directories) ) {
	if ( $dir ne "src" ) { $value++; }
        # Use a constant upward path, this allows arbitary source folder tree depth 
        $key = "../../../".$key.":";
	if ($add_vpath_level == 2) {
	    $key = "../../".$key;
	}
	push @vpath, $key ;
    }
    print MAKEFILE @vpath, "\n\n";

#__________________________________________________________________________________________________________________________________ 
# write compile and link information

    print MAKEFILE "%.o: %.cu\n",
        "\t\$(NVCC) \$(NVCFLAGS) -c \$<\n",
        "\n";

    print MAKEFILE "%.o: %.f\n",
        "\t\$(F77) \$(F77FLAGS) -c \$<\n",
        "\n";

    if (($dir =~ m/^externals/) && ($target =~ /^sx/)) {
	print MAKEFILE "\n",
            "%.o: %.f90\n",
	 "\t\$(FC) \$(FlibFLAGS) -c \$<\n",
	 "\n",
	 "%.o: %.F90\n",
	 "\t\$(FC) \$(FlibFLAGS) -c \$<\n";
    } else {	
	# Extra rule for JSBACH source files which need to be pre-processed by dsl4jsb.py
	print MAKEFILE "%_dsl4jsb.f90: %.f90\n",
	 "\t@ ../../../externals/jsbach/scripts/dsl4jsb/dsl4jsb.py -v -p _dsl4jsb -i \$<  -t .\n" ,
	 "\n",
	 "%.o: %.f90\n",
	 "\t\$(FC) \$(FFLAGS) -c \$<\n",
	 "\n",
	 "%.o: %.F90\n",
	 "\t\$(FC) \$(FFLAGS) -c \$<\n";
    }
    
#     print MAKEFILE "%.obj: %.f90\n";
#     print MAKEFILE "\t\$(FC) \$(FFLAGS) -c \$<\n";
#     print MAKEFILE "\n";
  
#__________________________________________________________________________________________________________________________________
# write all source files but not the program files
    
    my %seen = ();
    @seen{values(%target_programs)} = undef;
    
    my @sources = ();
    if ( "$dir" eq "support" ) {
        my ($cpu, $vendor, $os) = split /-/, $target;
        if ( "$os" eq "superux") {
            push (@sources, "rtc_sx.s");
        }

        #If we are linking against an external cdi lib, skip support/cdilib.c.
        if (exists($ENV{'CDIROOT'})
                && $ENV{'CDIROOT'} !~ '^\s*$') {
            $seen{"cdilib.c"} = undef;
        }
    }

    @sources = grep { $_ ne "version.c" and not exists($seen{$_}) } @source_files;

    my %unique = ();
    my @uniq_sources = grep { ! $unique{$_} ++ } @sources;    

    print MAKEFILE "SRCS =\t", IndentList(8, 0, \@uniq_sources),
     "\n\n";

    my @objects = ();
    foreach my $src (@uniq_sources) {
	$src =~ s/\.\S+$/\.o/;
	push @objects, $src;
    }
    print MAKEFILE "OBJS =\t", IndentList(8, 0, \@objects), "\n\n";

#__________________________________________________________________________________________________________________________________
# write targets

    my @target_all = ();
    while ( my ($key, $value) = each(%target_programs) ) {
	push @target_all, "../bin/$key";
    }
    if ( "$dir" ne "$iconrelsrcdir" ) {
	if ( $dir =~ m/^externals/) {
	    print MAKEFILE ".PHONY: \$(LIB)\n\n",
	     "all: \$(LIB)\n\n",
	     "\$(LIB): ../../../lib/lib\$(LIB).a\n";
	} else {
	    print MAKEFILE "all: \$(LIB)\n\n",
	     "\$(LIB): ../lib/lib\$(LIB).a\n";
	}
    } else {
	print MAKEFILE "all: create_version_c \$(OBJS) ",
            IndentList(13, 0, \@target_all);
    }
    print MAKEFILE "\n\n";

    if ( "$dir" ne "$iconrelsrcdir" ) {
	if ( $dir =~ m/^externals/) {
        print MAKEFILE "../../../lib/lib\$(LIB).a: \$(OBJS)\n",
            "\t\$(AR) \$(ARFLAGS) ../../../lib/lib\$(LIB).a \$(OBJS)\n",
            "\t\@for modfile in \$(wildcard *.mod); do \\\n",
            "\t\tcp \$\$modfile ../../../include; \\\n",
            "\t done\n\n";
        {
          if ( $dir =~ m/yac/) {
            print MAKEFILE 'CFLAGS += -I../../../../../', $dir,
                ' -I../../../../../', $dir, "/xml\n";
          }
          else {
            my $include_dir = $dir;
            $include_dir =~ s/src/include/;
            print MAKEFILE "CFLAGS += -I../../../../../$include_dir\n";
          }
        }
        print MAKEFILE "FFLAGS := \$(subst ../module,../../../module, \$(FFLAGS))\n";	    
        if ( $dir =~ m/self/) {
           print MAKEFILE 'FFLAGS := $(subst -C=all,,$(FFLAGS))';print MAKEFILE "\n";
        }
        print MAKEFILE "\n\n";
	} else {
	    print MAKEFILE "../lib/lib\$(LIB).a: \$(OBJS)\n",
	     "\t\$(AR) \$(ARFLAGS) ../lib/lib\$(LIB).a \$(OBJS)\n\n";
	}

	if ( "$dir" eq "support" ) {
        print MAKEFILE "ifeq (\$(ARCH), SX)\n",
         "rtc_sx.o: rtc_sx.s\n",
	     "\t\$(AS) -c rtc_sx.s\n",
	     "endif\n\n";
	}

    } else {
	print MAKEFILE <<"__EOF__"
.PHONY: create_version_c

create_version_c:
\t../../../config/pvcs.pl --srcdir ../../..

version.c: | create_version_c

version.o: version.c

libicon.a: \$(OBJS)
\t\$(AR) \$(ARFLAGS) \$@ \$(OBJS)

__EOF__
;
	while ( my ($key, $value) = each(%target_programs) ) {
	    my $okey = $key;
	    if ( "$okey" eq "jsb4_driver" ) {
	        $okey = "jsb4_driver_dsl4jsb.o" ;
	    } else {
	        $okey =~ s/ *$/.o/;
	    }
	    print MAKEFILE "$okey: $value
../bin/$key: $okey libicon.a version.o
\t\$(FC) \$(LDFLAGS) -o \$@ \$< libicon.a version.o \$(LIBS)

";
	}
    }
    
    print MAKEFILE "clean:\n",
     "\trm -f *.o *.mod ../module/*.mod\n",
     "\n";
    
#__________________________________________________________________________________________________________________________________
# print Fortran module dependencies

# don't need c implicit rules are used and don't like dependecies on system header files
#    @c_sources = grep /\.c$/, @sources;
#
#    if ( $#sources > 0 ) {
#        for my $file ( @c_sources ) {
#	    my ($object) = $file;
#	    $object =~ s/c$/o/;
#	    my (@includes) = @{$c_includes{$file}};
#	    print MAKEFILE "$object: $file ",
#	        IndentList (length($object)+length($file)+3, 0, @includes);
#	    print MAKEFILE "\n\n";
#	}
#    }
    
    if ( $dir =~ m/self/) {
      my @_myvpath = @vpath;
      shift @_myvpath;
      my $_myvpath = join('',@_myvpath);
      chop $_myvpath;
      print MAKEFILE "-include ",$_myvpath,"/../Makefile.depend";print MAKEFILE "\n";
    } else {
    for my $file ( keys %module_usage ) {
	next if $file =~ /c$/;
	my ($object) = $file;
	$object =~ s/f90$/o/;
        my %seen_module_usage = ();
	my (@modules) = grep { ! $seen_module_usage{$_} ++ } @{$module_usage{$file}};
	my (@dependencies) = ();
	for my $i ( 0 .. $#modules) {
	    my $ofile = "";
	    next if ! defined $module_definitions{$modules[$i]};
	    ($ofile) = $module_definitions{$modules[$i]};
	    next if ($ofile eq "");
	    $ofile =~ s/f90$/o/;
	    next if $object =~ $ofile;
	    push @dependencies, $ofile;
	}
	next if $object =~ $file;
	print MAKEFILE "$object: $file ",
            IndentList (length($object)+length($file)+3, 0, \@dependencies),
            "\n\n";
    } }

    close (MAKEFILE);
    
}
print "\n";
#__________________________________________________________________________________________________________________________________
#
exit;
#__________________________________________________________________________________________________________________________________
#
#__________________________________________________________________________________________________________________________________
#
# Subroutines:
#
#__________________________________________________________________________________________________________________________________
#
sub ScanDirectory {

# arguments to function    
    my ($workdir) = shift; 
    my ($workpath) = shift;
    my ($level)   = shift;    

# local work
    my($startdir) = &cwd;

    $level++;

    chdir($workdir) or die "Unable to enter dir $workdir:$!\n";

    my($pwd) = &cwd;
    $pwd =~ s/$prefix//;
    $pwd =~ s/^\///;

    $vpath_directories{$workpath} = $level;
    my %ignored_entries;
    {
      my @ignored_list = ('.', '..', '.svn', 'templates', 'hydro',
                          'interface', 'nh', 'phys', 'sw_options',
                          'ecrad');
      push(@ignored_list, 'ocean', 'sea_ice', 'hamocc')
          if ($enable_ocean eq "no");
      push(@ignored_list, 'lnd_phy_jsbach') if ($enable_jsbach eq "no");
      push(@ignored_list, 'serialization') if ($enable_serialization eq "no");
      push(@ignored_list, 'testbed') if (($enable_testbed eq "no")
                                             and ($workpath eq "$iconsrcdir"));
      @ignored_entries{@ignored_list} = undef;
    }
    opendir(DIR, ".") or die "Unable to open $workdir:$!\n";
    my @names = grep { !exists($ignored_entries{$_}) } readdir(DIR);
    closedir(DIR);

    my $suffix_list ="";

    foreach my $name (@names){

        if (-d $name){
            my $nextpath="$workpath/$name";
            ScanDirectory($name, $nextpath, $level);
        } else {
	    if ($enable_gpu eq "yes") {
		$suffix_list = '\.[c|f|F]{1}(u|90|95|03)?$';
	    } else {
		$suffix_list = '\.[c|f|F]{1}(90|95|03)?$';
	    }
	    if ($name =~ /$suffix_list/) {
                if ($workpath =~ "lnd_phy_jsbach") {
                    # For JSBACH, use the pre-processed source file located in the build directory
                    # These files need to be initially created by configure, with an additional "_dsl4jsb" before the suffix,
                    # so that the Makefile dependencies can be generated here.
                    my ($bname, $path, $suffix) = fileparse($name, '\.[^\.]*');  # parts of original source file
                    $name = $bname . "_dsl4jsb" . $suffix ;                      # name of pre-processed JSBACH files
                    open F, '<', $iconppdir . "/" . $name
                        or die("Cannot open file ${iconppdir}/${name}", $!);
                } else {
		    open F, '<', $name
                        or die("Cannot open file $name", $!);
                }
		my @lines = <F>;
		close (F);

		push @source_files, $name;

		my @filteredLines;
		simplifiedCPPFilter(\@lines, \@filteredLines);

		foreach (@filteredLines) {

		    if (/^ *MODULE/i && ! /procedure/i) {
			s/MODULE//i;
			s/\s//g;
			(my @d0) = split '!';
			$module_definitions{$d0[0]} = $name;
		    }
		    if ( (/^ *PROGRAM/i) && ($name !~ /\.c$/) ) {
			s/PROGRAM//i;
			s/\s//g;
			(my @d1) = split '!';
			$target_programs{$d1[0]} = $name;
		    }
		    if (/^ *USE/i) {
			s/USE//i;
                        (my @dt2) = split ('!', $_, 2);
			(my @d2) = split ',', $dt2[0];
			$d2[0] =~ s/\s//g;
			if (exists $module_usage{$name}) {
			    push @{ $module_usage{$name} }, $d2[0]; 
			} else {
			    $module_usage{$name}[0] = $d2[0];
			}
		    }
		    if (/^ *INCLUDE\s*[\'|\"]+/i) {
                        (my @dt3) = split ('!', $_, 2);
			my $d3 = $dt3[0];
			$d3 =~ s/\s*include\s*[\'|\"]+\s*(\w*)/$1/i;
			$d3 =~ s/(\'|\")//;
			$d3 =~ s/\s//g;
			if (exists $fortran_includes{$name}) {
			    push @{ $fortran_includes{$name} }, $d3; 
			} else {
			    $fortran_includes{$name}[0] = $d3;
			}
		    }
		    if (/^ *#include/) {
			if ( $name =~ /\.c$/ ) {
			    my $d4 = $_;
                            $d4 =~ s/\s*#include\s*(\w*)/$1/;
			    $d4 =~ s/\s//g;
			    if (exists $c_includes{$name}) {
				push @{ $c_includes{$name} }, $d4; 
			    } else {
				$c_includes{$name}[0] = $d4;
			    }
			} elsif ( $name =~ /\.[f|F]{1}(90|95|03)?$/ ) {
			    (my @d4) = split ' ';
			    $d4[1] =~ s/\s//g;
			    if (exists $fortran_includes{$name}) {
				push @{ $fortran_includes{$name} }, $d4[1]; 
			    } else {
				$fortran_includes{$name}[0] = $d4[1];
			    }
			}
		    }
		}
	    }
	}
    }
    chdir($startdir) or die "Unable to change to dir $startdir:$!\n";
}
#__________________________________________________________________________________________________________________________________
#
sub IndentList {

# arguments to function
    my ($columns) = 78 - shift(@_);
    my ($extratab) = shift(@_);
    my ($text) =  shift(@_);

    my ($wordlength);
    my @words;
# local work
    if (defined $$text[0]) {
      $columns -= length($$text[0]);
      foreach my $word (@$text) {
	$wordlength = length($word);
	if ($wordlength + 1 < $columns) {
          push(@words, " $word");
          $columns -= $wordlength + 1;
	}
	else {
          if ($extratab) {
            push(@words, " \\\n\t\t$word");
            $columns = 62 - $wordlength;
          }
          else {
            push(@words, " \\\n\t$word");
            $columns = 70 - $wordlength;
          }
	}
      }
    }
    return @words
}
#__________________________________________________________________________________________________________________________________
#
sub BuildSetup {
    
# arguments to function
    my ($build_path) = shift;
    my ($build_target) = shift;
    my ($build_directories) = shift;

# local work
    my @path7 = ();

    my ($path1) = "$build_path/build";
    my ($path2) = "$path1/$build_target";
    my ($path3) = "$path2/include";
    my ($path4) = "";
    if ( $build_target !~ /superux/ ) {
	$path4 = "$path2/module";
    }
    my ($path5) = "$path2/bin";
    my ($path6) = "$path2/lib";
    foreach my $dir ( @{$build_directories} ) {    
	push @path7, $path2."/".$dir;
    }

    if ( -d $build_path ) {
	if ( ! -d $path1 ) {
	    mkdir $path1 || die "Couldn't create build directory";
	}
	mkdir "$path2" || die "Couldn't create build directory";
    } else {
        die "Couldn't create build directory";
    }
    
    if ( ! -d $path3 ) {
	mkdir $path3 || die "Couldn't create module directory";		
    }
    if ( (! -d $path4) && ($path4 =~ /module/) ) {
	mkdir $path4 || die "Couldn't create bin directory";		
    }
    if ( ! -d $path5 ) {
	mkdir $path5 || die "Couldn't create lib directory";		
    }
    if ( ! -d $path6 ) {
	mkdir $path6 || die "Couldn't create lib directory";		
    }
    foreach my $path ( @path7 ) {        
	if ( ! -d $path ) {
	    mkdir $path || die "Couldn't create object directory";
	}		
    }

    return $path2;
}
#__________________________________________________________________________________________________________________________________
#
sub simplifiedCPPFilter {

    my ($lines, $filteredLines) = @_;

    #______________________________________________________________________________________________________________________________
    # regex for cpp #ifdef/#elif/#else/#endif (simple expressions only)

    my $ifdef_re = qr{^\s*#\s*ifdef\s+(\w*)};
    my $else_re = qr{^\s*#\s*else\b};
    my $elif_re = qr{^\s*#\s*elif\s+(\w*)};
    my $endif_re = qr{^\s*#\s*endif\b};

    #______________________________________________________________________________________________________________________________
    # input of a cpp definable variable not supposed to be part the source code

    my $exclude = "__COSMO__";
    my $printit = 1;           

    #______________________________________________________________________________________________________________________________
    # parse 

    if ( (not defined $exclude) || ($exclude =~ /^\s*$/) ) {
	$exclude = "__THIS_SHOULD_NEVER_EVER_HAPPEN__20140819__";
    }
    
    my $ifdef_count = 0;
    my @ifdef_stack;
   
    foreach (@$lines) {
	
	if ( /$ifdef_re/ ) {
	    $ifdef_count++; 
	    push @ifdef_stack, $1;
	    if ( /$exclude/ ) {
		push @$filteredLines, $_;
		$printit = 0;
	    }
	}
	
	if ( /$elif_re/ ) {
	    pop @ifdef_stack;
	    push @ifdef_stack, $1;
	    if ( /$exclude/ ) {
		push @$filteredLines, $_;
		$printit = 0;
	    } else {
		$printit = 1;            
	    }
	}
	
	if ( /$else_re/ ) {
	    if ( ($printit == 0) &&  ($ifdef_stack[-1] =~ m/$exclude/) ) {
		$printit = 1;
	    } 
	}
	
	if ( /$endif_re/ ) {
	    if ( ($printit == 0) && ($ifdef_stack[-1] =~ m/$exclude/) ) {
		$printit = 1;
	    }         
	    $ifdef_count--; 
	    pop @ifdef_stack;
	}
 
	push @$filteredLines, $_ if $printit;
    }
    return;
}
#__________________________________________________________________________________________________________________________________
