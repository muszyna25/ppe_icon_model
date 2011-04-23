#
# Parse an iconprotex conforming file and convert to doxygen documentation.
#
# Try to locate iconprotexs prolog sections,
# and try to move the contained information
# before the corresponding definitions.
# Besides, try to resolve LaTeX specific constructs,
# and try to guess formatting from typography
# (eg enumerations, brief and full description etc.)
#

use strict;

use File::Basename;
use POSIX qw(strftime);
use Getopt::Std;
$Getopt::Std::STANDARD_HELP_VERSION = 1;
$Getopt::Std::STANDARD_HELP_VERSION = $Getopt::Std::STANDARD_HELP_VERSION;

#
# Definition of constants
# 
 
#
# Subelements of prologue section.
#
my $prolog_key_subroutine_interface = "!SUBROUTINE INTERFACE:";
my $prolog_key_version_control = "!VERSION CONTROL:";
my $prolog_key_end = '--- End';
my @prolog_keys_ignored = (
    "!USES:",
    "!GLOBAL VARIABLES:",
    "!LOCAL VARIABLES:",
    "!INPUT PARAMETERS:",
    "!INPUT/OUTPUT PARAMETERS:",
    "!OUTPUT PARAMETERS:",
    "!RETURN VALUE:",
    "!INTERFACE:",
    "!FUNCTION INTERFACE:",
    "!PUBLIC TYPES:",
    "!PRIVATE TYPES:",
    "!PUBLIC MEMBER FUNCTIONS:",
    "!PRIVATE MEMBER FUNCTIONS:",
    "!PUBLIC DATA MEMBERS:",
    "!PARAMETERS:",
    "!ARGUMENTS:",
    "!DEFINED PARAMETERS:",
    );
my @prolog_keys = (
    "!REVISION HISTORY:",
    "!COPYRIGHT:",
    "!LICENSE:",
    "!WARRANTY:",
    "!BUGS:",
    "!SEE ALSO:",
    "!SYSTEM ROUTINES:",
    "!FILES USED:",
    "!REMARK:",
    "!REMARKS:",
    "!TO DO:",
    "!CALLING SEQUENCE:",
    "!AUTHOR:",
    "!CALLED FROM:",
    $prolog_key_version_control,
    $prolog_key_end,
    @prolog_keys_ignored,
);
# Re-format to yield regexp fragments.
my $prolog_pattern = join('|', @prolog_keys);
my $prolog_pattern_ignored = join('|', @prolog_keys_ignored);

#
# Declaration of program variables.
#

# Control options
my %command_line_options;
my $only_add_doxygen = 0;         # Set this to 0 if you want to delete
                                  # the old comments

# State variables
my $state = 'init';               # Scanner state
my @rescan_buffer = ();           # Buffer used for line re-scanning
                                  # so are independent of matching order

my @line_buffer = ();             # Buffer for postponed output

my $current_indent      = '';     # Indent for documentation comments,
                                  # taken from corresponding element.
my $current_title       = '';     # Main page title
my @current_authors     = ();     # Main page authors
my $current_date        = '';     # Main page date
my @current_description = ();     # Prologue description
my $current_prolog_key = '';      # Prologue subelement type 
my @current_prolog_elements = (); # Subelements as (type, text lines) pairs
my $current_prolog_type = '';     # Object type for this prologue
my $current_prolog_name = '';     # Object name for this prologue

my $line = '';

my $rescan_line = 0;

my $page_number = 0;

#
# Main routine
#

# Command line processing
getopts('a', \%command_line_options);
defined($command_line_options{a}) and $only_add_doxygen = 1;

while ( $line = @rescan_buffer ? shift(@rescan_buffer) : <> ) {

    if($only_add_doxygen) {
        # All lines are pushed into the line buffer,
        # unless they were scanned already.
        if(! $rescan_line) {
            push(@line_buffer, $line);
        }
    }

    #
    # State transitions
    #

    # By default, lines are just echoed.
    # When a program/module line is found, we start processing.
    if ( $state eq 'init' ) {
        if ( $line =~ /^s*!BOI/ ) {
            $state = 'intro';
        }
        elsif ( $line =~ /^s*!BOP/ ) {
            $state = 'prolog';
        }
        elsif ( $line =~ /^s*![BE]OC/ ) {
            
            # Ignored!

            # If we are in add only mode, we must avoid using the buffer
            # when in init state
            if($only_add_doxygen) {
                print($line);
                pop(@line_buffer);
            }
        }
        elsif ( $line =~ /^\s*($|[!#])/ ) {
            
            print($line);
            
            if($only_add_doxygen) {
                # This line was already printed,
                # so we do not need to buffer it.
                pop(@line_buffer);
            }
            
        }
        else {
            push( @rescan_buffer, $line );
            $state = 'program';
        }
    }

    # When we found a main program, buffer all lines
    # until we find documentation comments.
    elsif ( $state eq 'program' ) {
        if ( $line =~ /^\s*!BOI/ ) {
            $state = 'program/intro';
        }
        elsif ( $line =~ /^\s*!BOP/ ) {
            $state = 'prolog';
        }
        elsif ( $line =~ /^s*![BE]OC/ ) {
            
            # Ignored!
        }
        elsif ( is_eof($line) ) {

            # Flush pending output from file
            emit_program_text();
            $state = 'init';
        }
        else {

            # Save program or module name as default
            if($line =~ /(program|module)\s+(\w+)/i) {
                $current_prolog_type = $1;
                $current_prolog_name = $2;
            }

            ### DEBUG
            ### $line =~ s/((end\s+)?(program|module)\s+\w+)/${1}_debug/i;
            ### pop(@line_buffer);
            ### push( @line_buffer, $line );

            unless($only_add_doxygen) {
                push( @line_buffer, $line );
            }
        }
    }

    # When we found the introduction section,
    # we gather all infos for the title page.
    # There are two states, depending on where we must return to.
    elsif ( $state eq 'intro' ||
            $state eq 'program/intro' ) {
                
        if ( $line =~ /!TITLE:\s*(.*?)[\s.]*$/ ) {

            # Strip trailing whitespace and periods
            $current_title = $1;
        }
        elsif ( $line =~ /!AUTHORS:\s*(.*?)\s*$/ ) {
            push( @current_authors, [ $1, '' ] );
        }
        elsif ( $line =~ /!AFFILIATION:\s*(.*?)\s*$/ ) {
            $current_authors[$#current_authors]->[1] = $1;
        }
        elsif ( $line =~ /!DATE:\s*(.*?)\s*$/ ) {
            $current_date = $1;
        }
        elsif ( $line =~ /^\s*!EOI/ ) {
            emit_mainpage();
            if($state eq 'intro') {
                $state = 'init';
            }
            else {
                $state = 'program';
            }
        }
        else {

            # All other lines are echoed.
            unless($only_add_doxygen) {
                push( @line_buffer, $line );
            }
        }

    }

    # When in prologue, we need to process the text blocks
    elsif ( $state eq 'prolog' ) {

        $current_prolog_key = '';

        if ( $line =~ /$prolog_key_end|$prolog_pattern_ignored/o ) {
            # End marker may be ignored here.
            # In description and other states, end marker is significant.
        }
        elsif( $line =~ /$prolog_key_version_control/o ) {
            # We have to fake the version control element.
            # Formatting is done by emitter.
            push( @current_prolog_elements,
                  { key => $prolog_key_version_control,
                    text => [] } );
        }
        elsif( $line =~ /$prolog_key_subroutine_interface/o ) {
            $state = 'prolog/interface';
        }
        elsif ( $line =~ /!(MODULE|I?ROUTINE):\s*(\w+)/ ) {
            
            my ($type, $name) = ($1, $2);
            
            if($type eq 'MODULE' || $type eq 'ROUTINE') {
                
                $type eq 'MODULE' and $type = 'module';
                $type eq 'ROUTINE' and $type = 'program';
                
                if($current_prolog_type) {
                    if(lc($type) ne lc($current_prolog_type) ||
                       lc($name) ne lc($current_prolog_name)    )
                    {
                        warn("$ARGV: using $current_prolog_type $current_prolog_name instead of alleged $type $name\n");
                    }
                }
                else {
                    $current_prolog_type = $type;
                    $current_prolog_name = $name;
                }
            }
            else {
                $current_prolog_type = 'subroutine|function';
                $current_prolog_name = $name;
            }
        }
        elsif ( matches_rescan_blanked($line, '!DESCRIPTION:') ) {
            $state = 'prolog/description';
        }
        elsif (
            $current_prolog_key =
                matches_rescan_blanked($line, $prolog_pattern)
        ) {
            push( @current_prolog_elements, 
                  { key => $current_prolog_key,
                    text => [] } );
            $state = 'prolog/other'
        }
        elsif ( $line =~ /^\s*!([BE]O[PIC])/ ||
                is_eof($line)          ) {
            # We should have met an EOP, but...
            my $token = $1;
            if ( $1 ne 'EOP' ) {
                warn("$ARGV: prologue of $current_prolog_type $current_prolog_name is ended by $1 instead of EOP\n");
            }
            emit_prolog();
            $state = 'init';
        }
        elsif ( $line !~ /^\s*!/ ) {
            
            # Non-comment lines are echoed.
            unless($only_add_doxygen) {
                push( @line_buffer, $line );
            }
        } 
        else {

            # All other lines are echoed.
            unless($only_add_doxygen) {
                push( @line_buffer, $line );
            }
        }

    }

    # Also try to figure out signature from subroutine interface
    elsif ( $state eq 'prolog/interface' ) {
        if ( $line =~ /^\s*!([BE]O[PIC]|!*\s+(!DESCRIPTION:|$prolog_pattern|!))/o ||
             is_eof($line)                                                    ) {
            # Found next token; needs to be rescanned
            push( @rescan_buffer, $line );
            $state = 'prolog';
        }
        else {
            
            # If we have the signature, save values for anchoring
            if($line =~ /(subroutine|function)\s+(\w+)/i) {
                
                my ($type, $name) = ($1, $2);
                
                if($current_prolog_type) {
                    if($type !~ /($current_prolog_type)/i ||
                       lc($name) ne lc($current_prolog_name) )
                    {
                        warn("$ARGV: using $type $name instead of alleged $current_prolog_type $current_prolog_name\n");
                    }
                }
                $current_prolog_type = $type;
                $current_prolog_name = $name;
            }
            
            # All other lines are echoed.
            unless($only_add_doxygen) {
                push( @line_buffer, $line );
            }
        }
    }
    
    # Compile description or other elements in prolog section
    # and translate LaTeX specific syntax
    elsif ( $state eq 'prolog/description' ||
            $state eq 'prolog/other'          ) {

        if ( $line =~ /^\s*!([BE]O[PIC]|!*\s+(!DESCRIPTION:|$prolog_pattern|!))/o ||
             $line !~ /^\s*!/ ||
             $line =~ /^\s*$/ ||
             is_eof($line)                                                    ) {
            # Found next token; needs to be rescanned
            push( @rescan_buffer, $line );
            
            # Weed out LaTeX specific constructs.
            if( $state eq 'prolog/description') {
                delatexify(\@current_description);
            }
            else {
                delatexify($current_prolog_elements[$#current_prolog_elements]->{text});
            }
            
            $state = 'prolog';
        }
        elsif ( $line =~ /^\s*! ?(.*?)\s*$/ ) {
            $line = $1;
            if( $state eq 'prolog/description') {
                push( @current_description, $line );
            }
            else {
                my $buffer =
                    $current_prolog_elements[$#current_prolog_elements]->{text};
                push( @$buffer, $line );
            }
        }
        else {

            # All other lines are echoed.
            unless($only_add_doxygen) {
                push( @line_buffer, $line );
            }
        }
    }

    # Inject a special EOF token to avoid closing down 
    # everything outside of the loop.
    # While being passed down the state stack,
    # this line may be seen more than once, so make sure we only insert it once.
    if(eof && ! is_eof($line)) {
        push(@rescan_buffer, get_eof());
    }
    
    if($only_add_doxygen) {
        # For complete reproduction of the original file,
        # we need to know whether the next line was already scanned.
        $rescan_line = @rescan_buffer;
    }
    

}

#
# End of main routine
# 

#
# Subroutines
#

sub get_eof {
    '<<<<EOF>>>>';
}

sub is_eof {
    $_[0] eq get_eof();
}

sub emit_program_text {
    
    for $line (@line_buffer) {
        print $line;
    }
    
    # Empty flushed buffer
    @line_buffer = ();
}

sub emit_mainpage {
    my $line;

    print_begin();

    ++$page_number;
    my $file = lc(basename($ARGV, '.iconprotex'));
    $file =~ s/[^[:alnum:]]*//g;

    my $tag = "page$file$page_number";

    print_cont("\@page $tag $current_title");
    print_cont();

    # Author list
    for $line (@current_authors) {
        my ( $authors, $affiliation ) = @{$line};
        print_cont('@author');
        print_cont( '    ' . $authors );
        if ($affiliation) {
            print_cont( '    (' . $affiliation . ')' );
        }
        print_cont();
    }
    print_cont();
    
    # Date notice
    if($current_date eq '\today') {
        $current_date = strftime('%Y-%m-%d %H:%M:%S UTC', gmtime());
    }
    print_cont('@date '.$current_date);
    print_cont();
    
    print("\n");
}

sub emit_prolog {
    my $line;
    my $brief;
    my $anchor_found = 0;
    
    # Emit all buffered lines up to anchor line.
    while($line = shift(@line_buffer)) {
        if($line =~ /^(\s*).*?($current_prolog_type)\s+$current_prolog_name/i) {
            $current_indent = $1;
            unshift(@line_buffer, $line);
            $anchor_found = 1;
            last;
        }
        print($line);
    }
    
    # Safety catch...
    $anchor_found or
        die("Oops: missing anchor line ($current_prolog_type $current_prolog_name)\n");
        
    print_begin();
    
    emit_prolog_description();
        
    for my $element (@current_prolog_elements) {
        emit_prolog_element($element);
    }
    
    # Now print deferred output from file
    emit_program_text();

    # Clean up
    $current_prolog_type = '';
    $current_prolog_name = '';
    @current_prolog_elements = ();

}

sub emit_prolog_description {

    my $brief = get_brief_description();
    $brief .= '.' if $brief;
    
    emit_text_element('', $brief, \@current_description);

    @current_description = ();
}

sub emit_prolog_element {
    
    my $key = $_[0]->{key};
    
    ### # Weed out unwanted documentation elements.
    ### if($key =~ /$prolog_pattern_ignored/o) {
    ###     return;
    ### }
    
    # Special handling for special elements.
    if($key eq $prolog_key_version_control) {
        print_cont('$Id: n/a$');
        print_cont();
        return;
    }
    
    my $title = format_key($key);
    my $text = $_[0]->{text};

    emit_text_element($title, '', $text);
}

sub emit_text_element {
    my $title = $_[0];
    my $brief = $_[1];
    my $text  = $_[2];
    
    my $need_paragraph = 0;
    my $blank_lines_found = 0;
    my $enumeration_found = 0;
    
    if($title) {
        print_cont('@par '.$title);
    }
    if($brief) {
        print_cont($brief);
        print_cont($title ? '@par' : '');    
    }
    for my $line (@$text) {
        if($line =~ /^\s*$/) {
            $blank_lines_found = 1;
            if($enumeration_found) {
                print_cont('</ol>');
                $enumeration_found = 0;
            }
        }
        elsif($line =~ /^\s*\d+\.(\s.*)$/) {
            if(! $enumeration_found) {
              $enumeration_found = 1;
              print_cont('<ol>');
            }
            print_cont('<li>', $1)
        }
        elsif($enumeration_found && $line =~ /^\s+/) {
            print_cont($line);
        }
        else {
            if($enumeration_found) {
                print_cont('</ol>');
                $enumeration_found = 0;
            }
            if($blank_lines_found) {
                if($need_paragraph) {
                    print_cont($title ? '@par' : '');
                }
                else {
                    $need_paragraph = 1;
                }
                $blank_lines_found = 0;
            }
            print_cont($line);
        }
    }
    print_cont();
}

#
# Extract brief description from current desription text.
# For this, get text up to first period or any full line
# until the brief description has reached defined max length.
#
# (!) SIDE EFFECT: changes @current_description
#
sub get_brief_description {
    
    my $max_length = 76;
    
    my $todo = 1;   # State flag
    my $line = '';  # Current input line
    my $brief = ''; # Result variable
    my $sep = '';   # Separator for concatenation

    # Dissect brief description from full text.
    while($todo && @current_description) {
        
        $line = shift(@current_description);
        
        # Check for nasty combinations, like LaTeX environments.
        # In this case, split line and insist on finishing that line.
        my ($first, $terminator, $remainder) =
            split(/(\.\s*|\@f[{\[])/, $line, 2);
        defined($terminator) or $terminator = '';
        defined($remainder) or $remainder = '';
        if($terminator) {
            if(length($brief)+length($sep)+length($first) > $max_length) {
                # Don't use current line if brief overflowed.
                # But to make sure that the first sentence is OK,
                # we must repeat brief description!
                unshift(@current_description, $line);
                # Just in case that even the first line was too long...
                if(! $brief) {
                    $brief = $first;
                }
                else {
                    unshift(@current_description, $brief);
                }
            }
            else {
                # Otherwise, extend to use line up to period.
                $brief .= $sep.$first;
                if(substr($terminator, 0, 1) eq '.') {
                    # Leave out period and subsequent white space
                    unshift(@current_description, $remainder);
                }
                else {
                    # LaTeX may be inside a sentence,
                    # so - like above - we must repeat brief description!
                    unshift(@current_description, $terminator.$remainder);
                    unshift(@current_description, $brief);
                }
            }
            # In any case, we found a period, so stop here!
            $todo = 0;
        }
        else {
            if(length($brief)+length($sep)+length($line) > $max_length) {
                # Don't use current line if brief overflowed.
                # But to make sure that the first sentence is OK,
                # we must repeat brief description!
                # Besides, we now may stop here.
                unshift(@current_description, $line);
                # Just in case that even the first line was too long...
                if(! $brief) {
                    $brief = $line;
                }
                else {
                    unshift(@current_description, $brief);
                }
                $todo = 0;
            }
            else {
                # Here, we haven't found a period,
                # and brief would not overflow, so let's go on!
                $line and $brief .= $sep.$line;
                $brief and $sep = ' ';
            }
        }
    }

    my @tokens = ($brief =~ /\@f\$/g);
    if(@tokens % 2) {
        $brief .= '@f$';
    }
    
    return $brief;
}

# Emit the first line of a pre-positioned documentation block
sub print_begin {
    print $current_indent, '!>', "\n";
}

# Emit a continuation line for a pre-positioned documentation block 
sub print_cont {
    print $current_indent, '!! ', @_, "\n";
}

# Try to match a line with a given keyword.
# If the line matches, blank out the key and provide the result for rescanning.
sub matches_rescan_blanked {
    if($_[0] =~ /($_[1])/) {
        my $found = $1;
        my $blank = ' ' x length($found);
        $_[0] =~ s/$found/$blank/;
        push( @rescan_buffer, $_[0] );
        return $found;
    }
    return "";
}

# Re-format the ProTeX keywords
# from all capital letters to initial capitalization.
# Also remove leading and trailing punctuation.
sub format_key {
    my $key = shift(@_);
    $key =~ s/^!|:$//g;
    $key =~ s/(\w+)/ucfirst(lc($1))/ge;
    return $key;
}

# Try to replace LaTeX specific constructs by doxygen's equivalents.
sub delatexify {
    my $text = $_[0];
    
    my $line = '';
    my $new_line = '';
    my $token = '';
    my $payload = '';

    my $latex_mode = '';
    
    for my $i (0..$#$text) {
        $line = $text->[$i];
        $new_line = '';
        while($line =~ /\G(.*?)(\$\$?|\\item|\\(begin|end){.*?})/gc) {
            
            $payload = $1;
            $token = $2;
    
            #
            # Process payload in dependence on mode.
            #
            map_text($latex_mode, $payload);
            $new_line .= $payload;
            
            #
            # Process tokens
            #
            if(! $latex_mode) {
                if($token eq '$') {
                    $new_line .= '@f$';
                    $latex_mode = 'math';
                }
                elsif($token eq '$$') {
                    $new_line .= '@f[';
                    $latex_mode = 'equation';
                }
                elsif($token eq '\item') {
                    $new_line .= '<li> ';
                }
                elsif($token eq '\begin{itemize}') {
                    $new_line .= '<ul>';
                }
                elsif($token eq '\end{itemize}') {
                    $new_line .= '</ul>';
                }
                elsif($token eq '\begin{enumerate}') {
                    $new_line .= '<ol>';
                }
                elsif($token eq '\end{enumerate}') {
                    $new_line .= '</ol>';
                }
                elsif($token eq '\begin{verbatim}') {
                    $new_line .= '@verbatim';
                    $latex_mode = 'verbatim';
                }
                else {
                    $token =~ /(begin|end){(.*)}/;
                    if($1 eq 'begin') {
                        $latex_mode = $2;
                        $new_line .= '@f{'.$latex_mode.'}{';
                    }
                    else {
                        $new_line .= $token;
                    }
                }
            }
            else {
                if($latex_mode eq 'math' && $token eq '$') {
                    $new_line .= '@f$';
                    $latex_mode = '';
                }
                elsif($latex_mode eq 'equation' && $token eq '$$') {
                    $new_line .= '@f]';
                    $latex_mode = '';
                }
                elsif($latex_mode eq 'verbatim' && $token eq '\end{verbatim}') {
                    $new_line .= '@endverbatim';
                    $latex_mode = '';
                }
                elsif($token eq '\end{'.$latex_mode.'}') {
                    $new_line .= '@f}';
                    $latex_mode = '';
                }
                else {
                    $new_line .= $token;
                }
            }
        }
        if(defined(pos($line))) {
            $payload = substr($line, pos($line));
        }
        else {
            $payload = $line;
        }
        map_text($latex_mode, $payload);
        $new_line .= $payload;
        $text->[$i] = $new_line;
    }
    
}

sub map_text {
    my $latex_mode = shift(@_);
    if($latex_mode) {
        unless($latex_mode eq 'verbatim') {
            map_latex($_[0]);
        }
    }
    else {
        map_doxygen($_[0]);
    }
}

sub map_latex {
    # There seems to be a doxygen bug with end-of-line \\ instructions.
    # This is a workaround
    $_[0] =~ s/\\\\$/\\\\\\/;
}

sub map_doxygen {
    # Replace font formatting
    $_[0] =~ s/\\textit{(.*?)}/<i>$1<\/i>/g;
    $_[0] =~ s/{\s*\\it\s+(.*?)}/<i>$1<\/i>/g;
    $_[0] =~ s/{\s*\\bf\s+(.*?)}/<b>$1<\/b>/g;
    # Replace non-breaking spaces
    $_[0] =~ s/~/&nbsp;/g;
    # Backslashes embedded in identifiers are supposed to be LaTeX workarounds
    $_[0] =~ s/(\w)\\(\w)/$1$2/g;
    $_[0] =~ s/\\\\(\[.*?\])?/<br>/g;
    # Any remaining backslashes are quoted. Just to be sure.
    $_[0] =~ s/\\/\\\\/g;
}

