<!--
This file is written using Markdown language, which might make it difficult to
read it in a plain text editor. Please, visit ICON project page on DKRZ GitLab
(https://gitlab.dkrz.de/icon/icon/-/tree/master/config/buildbot/dwd_nec_patches)
to see this file rendered or use a Markdown viewer of your choice
(https://www.google.com/search?q=markdown+viewer).
-->

This directory contains a patch that we have (hopefully, temporarily) to apply
to YAXT on the NEC machine at DWD:

1. [`yaxt_mpi_abort_test.patch`](yaxt_mpi_abort_test.patch) &mdash; modifies the
`test_xmap_all2all_fail_run` test of YAXT, which expects a failure from the
`mpirun` process with exit code `3`. The existing workaround implemented in YAXT
for such cases is unfortunately not enough. Therefore, the test is modified to
pass if the `mpirun` program terminates with any non-zero exit code. See
https://gitlab.dkrz.de/icon/icon-cimd/-/merge_requests/36#note_72444 for more
information.
