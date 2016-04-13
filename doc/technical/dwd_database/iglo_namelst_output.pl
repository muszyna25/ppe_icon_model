#!/usr/bin/perl -w

=head1 NAME

B<iglo_namelst_output.pl> - Write model level variables to iglo namelists

=head1 SYNOPSIS

B<iglo_namelst_output.pl> [B<-d> F<dir>] [B<-D>]  F<nml_file1>]  [F<nml_file2> ...]
                 [B<-h>]

=head1 DESCRIPTION

B<iglo_namelst_output.pl> reads the tex-lists from the ICON Database
Reference Manual and sets ml_varlist in namelist files F<nml_file.>
for the C<iglo> script.

=head1 OPTIONS

=over 4

=item B<-d> F<dir>

Directory for the new namelist files

=item B<-h>

This help text.

=item B<-D>

Debug options. Writes extra output to STDERR.

=back

=head2 created by Helmut P. Frank, DWD, FE 13, 25.06.2015

=head2 Note

=head2 Modifications

26.11.2015

Adapt for tiles output and SNOWC.

=cut

use open qw< :encoding(ASCII) >;
use Cwd qw(abs_path);

my @in_files;
my $dir = '.';
my $debug;

while ( $arg  = shift @ARGV) {
    exec "perldoc $0" if ( $arg eq '-h');
    if ( $arg eq '-d') {
        $dir = shift @ARGV;
        next;
    }
    if ( $arg eq '-D') {
        $debug++;
        next;
    }
    if ( -r $arg) {
        push  @in_files, $arg;
    } else {
        print STDERR "File $arg is not readable!\n";
    }
}

if ( scalar( @in_files)==0 ) {
    print STDERR "No namelist files given!\n";
    exec "perldoc $0";
}

if ( ! -d $dir) {
    mkdir $dir or die "\nError! Could not create output directory $dir\n\n";
} elsif ( ! -w $dir) {
    die "\nError! Output directory $dir not writeable by you!\n\n";
}

my @tex_lists = ( "output_GRIB2_vv0.tex",
                  "output_GRIB2_forecast.tex",
                  "output_GRIB2_forecast_singlelvl.tex",
                  "output_GRIB2_nest.tex",
                  "analysis.tex");
my ( $vv, $grid, $nest, $ext);
my ( $var, @vars, $nml_file);
my ( %nml_var);
my %grid_vars;
my @out;
# Read the lists for the different grids
for $tex_list ( @tex_lists) {
    print STDERR "Read $tex_list\n";
    if ( $tex_list =~ /analysis/) {
        my @cols;
        open ANA, "<$tex_list";
        $out[0] = "list of variables for iau\n";
        $out[1] = '';
        while ( $_ = <ANA>) {
            @cols = split( /&/, $_);
            if ( scalar(@cols) > 4 && $cols[2] =~ /\bFG\b|\bAN\b/) {
                $cols[0] =~ s/ *//g;
                $cols[0] =~ s/\\_/_/g;
                $cols[0] =~ s/,/ /g;
                $cols[0] =~ s/\(\S+\)//;
                $cols[0] =~ s/\\footnote.*//;
                if ( $cols[-1] =~ /tiles:\s*$/) {
                    $tiles = 'tiles:';
                } else {
                    $tiles = '';
                }
                foreach $var ( split( ' ', $cols[0]) ) {
                    $out[1] .= "'$tiles$var', ";
                }
            }
        }
        close ANA;
#       The SMA needs 3h first guess data from the assimilation of the
#       following variables
        for $var ( qw(T_2M TD_2M U_10M V_10M) ) {
            $out[1] .= "'$var', "  unless ( $out[1] =~ /\b$var\b/);
        }
    } else {
        @out = qx(python build_varlists.py $tex_list);
        if ( $tex_list =~ /vv0.tex/) {
            $vv = '_0';
        } else {
            $vv = ''
        }
    }
    while ($_ = shift @out) {
        print STDERR  if ($debug);
        next if (/^\s*$/);
        if (/list of variables for/) {
            $grid = '';
            $nest = '';
            $ext  = '';

            $grid = '.main' if ( /for triangular/);
            $grid = '_ll'   if ( /for lon-lat/);
            $nest = ''      if ( /global grid/);
            $nest = '_eu'   if ( /local grid/);
            $ext  = '.i2l'  if ( /for int2lm/);
            $ext  = '.oth'  if ( /for triangular/);
            $ext  = '.iau'  if ( /for iau/);
#           Read the variables
            $_ = shift @out;
            print STDERR  if ($debug);
            s/[\[\]\s]//g;
            @vars = split( /,/, $_);
            $nml_file = 'namelst.output' . $nest . $grid . $vv . $ext;
            $nml_file =~ s/\.main_0.oth/_0/;
            $nml_file =~ s/\.main_ll.oth/_ll/;
            $nml_file = 'namelst.output.iau'  if ( $ext eq '.iau');
            print STDERR "\tnml_file=$nml_file\n";
            for $var (@vars) {
#               coordinates are written via output_grid=.true.
#               Do not include them in ml_varlist
                if ( $var =~ /\b[CEVR]LON\b|\b[CEVR]LAT\b/) {
                    $grid_vars{$var}++;
                    next;
                }

#               10 m winds on regular grid are written in special namelists
#               next if ( $grid eq '_ll' && $var =~ /[UV].*_10M\b/);
                if ( $grid eq '_ll' && $var =~ /[UV].*_10M\b/) {
                    $ml_var{"$nml_file.wind_10m"}{$var}++;
                    next;
                }
#               10 m gusts on global tringular grid are written hourly in special namelist
                if ( $grid ne '_ll' && $nest eq '' && $ext eq '.oth' && $var =~ /VMAX_10M\b/) {
                    $ml_var{"$nml_file.wind_10m"}{$var}++;
                    next;
                }

                $ml_var{$nml_file}{$var}++;
            }
        }
    }
}
# Special for grid coordinates
for $var ( qw(CLON CLAT ELON ELAT VLON VLAT) ) {
    $vv = "'$var'";
    unless ( defined($grid_vars{$vv}) && $grid_vars{$vv} > 0) {
        $ml_var{"namelst.output_0"}{"'-grid:$var'"} = 1;
        $ml_var{"namelst.output_eu_0"}{"'-grid:$var'"} = 1;
    }
}

# write the namelists
my $in  = undef;
my $out = undef;
my ( $in_file, $out_file, $nml_id);
my ( $in_abs, $out_abs);
my ( @ml_varlist, $n_varlist, $l);
my $output_start;
for $in_file (@in_files) {
    ( $nml_file = $in_file ) =~ s{.*/}{};

#   define name of output file.
    $in_abs   = abs_path( $in_file);
    $out_file = "$dir/$nml_file";
    $out_abs  = abs_path( $out_file);
    die "Error! Change output directory $dir!\n" .
        "Output file $out_file would overwrite\n" .
        " input file $in_file !\n"  if ( $out_abs eq $in_abs);

    if ( -l $in_file) {
        my $ilink = readlink $in_file;
        $ilink =~ s{.*/}{};
        print STDERR "Link $ilink -> $out_file\n";
#       symlink "$dir/$ilink", $out_file;
        symlink "$ilink", $out_file;
        next;
    }

    open( $in, "<", $in_file) or die "Error opening file $in_file!\n";
    print STDERR "Read $in_file\n" if ($debug);

    open( $out, ">", $out_file) or die "Error opening file $out_file!\n";
    print STDERR "Write $out_file\n";

    $nml_id = $nml_file;
    $nml_id =~ s/(namelst.output.*)(\.ass\.)/$1.main./;
    $nml_id =~ s/(namelst.output.*)(\.pre\.)/$1.main./;
    $nml_id =~ s/_eu// if ( $nml_id =~ /\.i2l$|\.iau$/);
    @ml_varlist = make_ml_varlist( $nml_id);

#   Exclude WW, DTKE_CON, ... from assimilation namelists
    if ( $nml_file =~ /\.ass\.|\.pre\./) {
        my @no_ass = grep( !/\bWW\b|\bDTKE_CON\b|\bH\w+_CON\b|\bHTOP_DC\b|\bHZEROCL\b|\bTC[HM]\b|\bASWDI|\bTQ._DIA\b/,
                           @ml_varlist);
        @ml_varlist = @no_ass;
    }

    $n_varlist = scalar( @ml_varlist);
    print STDERR "$nml_id: ml_varlist = @ml_varlist\n" if ($debug);

    while (<$in>) {
        if ( /output_bounds\s*=\s*(\S+)/) {
            $output_start = $1;
            if ( $output_start =~ /^(\d+)/) {
                $output_start = $1;
            } else {
                $output_start = 0;
            }
        }
#       Replace ml_varlist by the new list
        if ( $n_varlist && /^\s*ml_varlist *=/) {
          do { $_ = <$in>;
             } until ( /[^!]*\w+\s*=/);
             $l = " ml_varlist           = ";
             foreach $v (@ml_varlist) {
                 $var = $v;
                 if ( $output_start == 5400) {
#                    P, T, U, V are not necessary in the first guest output after 5400 s = 90 min
#                    because they are not prognostic variables
                     next if ( $var =~ /\bP\b|\bT\b|\bU\b|\bV\b/);
                 } else {
                     next if ( $var =~ /SNOWC/);
#  presently, tiles and SNOWC are only written to the first guess file forecast for 5400 s
                     $var =~ s/tiles://;
                 }
                 $l .= "$var, ";
                 if ( length($l) > 72) {
                     print $out "$l\n";
                     $l = ' ' x 24;
                 }
             }
             print $out "$l\n" if ( $l =~ /\w/);
             print { $out} $_;
        } else {
            print { $out} $_;
        }
    }
    close $in;
    close $out;
}
exit;

sub make_ml_varlist {
    my ( $nfile, $debug) = @_;
    my $r_var = $ml_var{$nfile};
    my @nml_vars = sort ignoreTiles keys %$r_var;
    return @nml_vars;
}

sub ignoreTiles {
    my ( $aa, $bb);
    ( $aa = $a ) =~ s/^'tiles:/'/;
    ( $bb = $b ) =~ s/^'tiles:/'/;
    $aa cmp $bb;
}
