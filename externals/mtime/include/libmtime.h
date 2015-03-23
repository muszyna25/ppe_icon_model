/** @file libmtime.h Documentation header file */

/** @mainpage libmtime - A model time management library
 *
 * This project aims to provide time, calenendar, and event handling
 * for model applications.
 *
 * libmtime provides some simple and well documented data structures of
 * common use:
 *
 * @arg calendar    the base calendar, either proleptic Gregorian, 365 day years, or 360 day years;
 * @arg date        a calendars date;
 * @arg time        a time at a day;
 * @arg datetime    date and time combined;
 * @arg timedelta   a difference in datetime;
 * @arg julianday   the Julian day to define the base time axis 
 *
 * @section CompileInstall Compilation and installation
 *
 * You can obtain the source code tarball on site
 *
 * This software uses the GNU Autotools build system: in order to
 * compile and install the application the standard installation
 * method must be used:
 @verbatim
 ./configure
 make
 make install
 @endverbatim
 *
 * More information is provided by the INSTALL file, which is
 * included in the software package.
 *
 * @section Usage
 *
 * The doxygen documentation you are reading is pretty rich. You
 * however can take a glance at the libmtime/test directory, which
 * contains some usage examples.
 *
 * @section License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section Dependency
 *
 * The library should work correctly without any dependency under any
 * Posix-compliant operating system. If you experience some issue
 * please contact me. The AUTHORS file provided with the software
 * distribution contain developers contacts.
 *
 */

/** @defgroup libmtime
 *
 * This library provides time, calenendar, and event handling
 * for model applications.
 */

/** @defgroup CBindings C language bindings
* @ingroup libmtime
*
* This module documents the C language bindings
*/

/** @defgroup FortranBindings Fortran language bindings
* @ingroup libmtime
*
* This module documents the Fortran language bindings
*/

