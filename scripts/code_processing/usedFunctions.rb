#!/usr/bin/env ruby

iconSourceDir, executable, gmonOutput, rest = ARGV

if [iconSourceDir,gmonOutput,executable].include?(nil) then
  warn "Cannot read iconSourceDir of gmonOutput file"
  warn "use #{__FILE__} iconSourceDir executable gmonOutput <ACTION=edit>"
  exit(1)
end
# =============================================================================
# SCAN {{{
# scan relevant source code
fileList       = Dir.glob("#{iconSourceDir}/**/*.f90")
scanSource     = IO.popen("ctags -x  --fortran-kinds=fs ../../src/ocean/*/*f90").readlines.map {|l| l.chomp.gsub(/\s +/,' ').split(' ')}
filesOfFunx    = {}; scanSource.each {|l| filesOfFunx[l[0]]    = l[3]}
patternsOfFunx = {}; scanSource.each {|l| patternsOfFunx[l[0]] = l[4..-1].join(' ')}
allFunx        = filesOfFunx.keys.sort.uniq

# scan the binary
calledFunx     = IO.popen("gprof -b -p #{executable} #{gmonOutput}").readlines.map(&:chomp).map {|l| l.split('MOD_')[-1] if /_MOD/.match(l)}.delete_if {|v| v.nil?}

# compute the intersection: only routines from the relevant source code are shown
usedFunx       = allFunx & calledFunx
puts usedFunx.join("\n")
#puts allFunx
#puts calledFunx.join("\n")
# }}}
# =============================================================================
# EDIT FILES {{{
options = {}
rest.each {|opt| k,v = opt.split('='); options[k] = v}
if options.has_key?('ACTION') and options['ACTION'] == 'edit' then
  usedFunx.each {|func|
    # find out the source code file where the routine is defined
    file    = filesOfFunx[func]
    pattern = patternsOfFunx[func]
    # Add <OptimizeEssential> in the line before
    cmd = "sed -i '/^ \\+#{pattern}/i   !<OptimizeEssential>' #{file}"
    puts cmd
    system(cmd)
  }
end
# }}}
#
# vim:fdm=marker

