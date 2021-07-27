# Release notes for icon-2.6.3

The new release of icon is available.

It consists of

- a large number of bugfixes, refactorings, and optimizations for the models infrastructure
- bugfixes and improvements in the build environment
- consolidation of all QUBICC based enhancements and bugfixes (for the sapphire physics)
- many add-ons for the sapphire physics port to GPUs based on OpenACC and many further steps on porting the NWP physics to GPUs
- ART has made its first step to an git submodule external (draft implementation not to be used yet- no warranty)
- fixes and improvements in the ocean code including hamocc
- improvements in the data assimilation NWP physics coupling
- tuning of data assimilation 
- added rrtm-gp as radiation scheme for the sapphire physics on GPU
- much progress on cdi-pio use in many icon components
- refactoring of mpi communication library (with a focus on GPU to GPU communication)

For many more details visit:

https://gitlab.dkrz.de/icon/wiki/-/wikis/Protocol-of-Release-Commits

June 7th, 2021

