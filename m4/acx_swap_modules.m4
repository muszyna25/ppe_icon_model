dnl acx_swap_modules.m4 --- swap and load user added defined modules 
dnl
dnl Copyright  (C)  2017  Luis Kornblueh <luis.kornblueh@mpimet.mpg.de>
dnl
dnl Version: 1
dnl Keywords:
dnl Author: Luis Kornblueh <luis.kornblueh@mpimet.mpg.de>
dnl Maintainer: Luis Kornblueh <luis.kornblueh@mpimet.mpg.de>
dnl URL: https://code.zmaw.de/projects/icon
dnl
dnl Redistribution and use in source and binary forms, with or without
dnl modification, are  permitted provided that the following conditions are
dnl met:
dnl
dnl Redistributions of source code must retain the above copyright notice,
dnl this list of conditions and the following disclaimer.
dnl
dnl Redistributions in binary form must reproduce the above copyright
dnl notice, this list of conditions and the following disclaimer in the
dnl documentation and/or other materials provided with the distribution.
dnl
dnl Neither the name of the Max PLanck Institute for Meteorology nor the
dnl names of its contributors may be used to endorse or promote products
dnl derived from this software without specific prior written permission.
dnl
dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
dnl IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
dnl TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
dnl PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
dnl OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
dnl EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
dnl PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
dnl PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
dnl LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
dnl NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
dnl SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
dnl
dnl
AC_DEFUN([ACX_SWAP_MODULES],
         [AC_REQUIRE([AC_PROG_AWK])
          AC_REQUIRE([AC_PROG_GREP])
          AC_REQUIRE([AC_PROG_SED])

          moduleshome_in_use="${MODULESHOME-}"
          load_modules_in_use=m4_default([$1],[""])
          shell_in_use="${SHELL##*/}"

          # in case specific modules are required do replace those only
                 
          AS_IF([test "x$load_modules_in_use" != "x"],   
                [AS_IF([test ! -z $moduleshome_in_use],
                       [$GREP TCLSH $moduleshome_in_use/init/$shell_in_use > /dev/null && { export TCLSH=$(which tclsh); } 
                        # load functions for module (C/tcl versions supported)
                        . $moduleshome_in_use/init/$shell_in_use],
                       [AC_MSG_ERROR([modules required, but module command not available!], [1])])

                 # list of modules loaded      
                 loaded_module_list=$(module list 2>&1 | $GREP -v Currently)
                 AC_MSG_NOTICE([Modules loaded from environment])
                 AS_ECHO(["$loaded_module_list"])])
		 
		 # swap PrgEnv if on a Cray
                 prgenv_loaded=$(module list -t 2>&1 | $AWK '/PrgEnv/{print}')	
		 prgenv_to_be_loaded=$(echo $load_modules_in_use | $GREP -o '\bPrgEnv[-\w]*')
                 AS_IF([test ! -z "$prgenv_to_be_loaded"],
                       [module swap $prgenv_loaded $prgenv_to_be_loaded]) 

		 modules_loaded=$(module list -t 2>&1 | $AWK '!/PrgEnv/{print}')
                 modules_to_be_loaded=$(echo $load_modules_in_use | $SED 's/\bPrgEnv[\S]*//')

		 for replacement in $modules_to_be_loaded
                 do		 
                     is_loaded=false 
                     for mod in $modules_loaded
                     do
		         AS_IF([test ${mod%%/*} = ${replacement%%/*}],
 		               [AS_IF([test "$mod" != "$replacement"],
                                      [AS_IF([test ! -z prgenv_loaded],
                                             [module unload $mod
                                              module load $replacement],
                                             [module swap $mod $replacement])],
                                      [is_loaded=true])])

		     done
                     AS_IF([test $is_loaded = false],
                           [module load $replacement])
                 done

                 # list of modules loaded finally for compiling     
                 loaded_module_list=$(module list 2>&1 | $GREP -v Currently)
                 AC_MSG_NOTICE([Modules loaded for compiling])
                 AS_ECHO(["$loaded_module_list"])])
])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://code.zmaw.de/projects/icon/wiki/ICON_licensees"
dnl license-default: "bsd"
dnl End:
