------------------------------------------------------
DWD Setup for BuildBot
------------------------------------------------------

[2017-07-31] : F. Prill, DWD

------------------------------------------------------
    Running the DWD Tests
    Creating New DWD Tests
        Nomenclature
        Queue Restrictions for User "buildbot"
        Registering Tests in the "dwd" Experiment List
        Test Data and Output Directory
        Signalling "SUCCESS" or "FAILURE"
    Admin Information
        User "buildbot"
        Access Rights
        Directory Structure for BuildBot at DWD
------------------------------------------------------




h1. DWD Setup for BuildBot

The "dwd" experiment list contains a set of medium-size test cases which

* are executed on the DWD machine only (all other builders inactive),
* consist of plain PBSpro run scripts and do not use the usual BuildBot macros.


h2. Running the DWD Tests

Select the following parameters on the Buildbot web site (https://buildbot.zmaw.de/icon/builders):

| Repository: | git@git.mpimet.mpg.de:icon-nwp.git |
| Properties: | Name: "list" ; value: "dwd"        |


h2. Creating New DWD Tests

The DWD experiment list "dwd" uses the "plain run script mechanism" that has been implemented by Leonidas and Ralf.
When importing existing scripts, only a few adjustments are necessary.


h3. Nomenclature

Run scripts are located in <pre>run/checksuite.xce.dwd.de</pre>

# It is necessary to rename the test scripts according to the following scheme: <pre>exp.<test name>.run</pre>
# It is necessary to set the stdout log output in the PBS header according to the following scheme: 
<pre>
#PBS -j oe
#PBS -o LOG.exp.<test name>.run.run.o 
</pre>

The scripts are started from the subdirectory @run@ with the help of a soft link to the @run/checksuite.xce.dwd.de@ subdirectory.
Therefore the environment variable @${PBS_O_WORKDIR}@ points to the @run@ subdirectory. 
This must be taken into account when using @${PBS_O_WORKDIR}@ together with relative paths.


h3. Queue Restrictions for User "buildbot"

Currently, buildbot tests on the DWD site are limited to 

* 4 concurrent runs for the queues @xc_norm_h@ and @xc_norm_b@
* 600 CPUs for the queue @xc_norm_h@ and 900 CPUs for the queue @xc_norm_b@.
  This corresponds to 25 nodes.

The current settings can be viewed by typing
<pre>
qstat -Q -f <queue name>
</pre>


h3. Registering Tests in the "dwd" Experiment List

New tests can be registered by typing
<pre>
cd scripts/buildbot_scripts/
addexp --list dwd checksuite.xce.dwd.de/<test name>
</pre>
Note that neither the prefix "exp." nor the suffix ".run" are used here!

This modifies the file @scripts/buildbot_scripts/experiment_lists/dwd@ which then has to be committed to the DKRZ repository.


h3. Test Data and Output Directory

Test data (initial state, boundary data etc.) may be loaded from the directory <pre>/lustre2/rwork0/routfor/test/icon/Checksuite_data</pre>.
Please ask T. Hanisch, e.g., for setting up a new subdirectory for your test case.

Output directories for the experiments may be created in @/lustre2/gtmp/buildbot@ (available via the environment variable @$TMPDIR@) where they are automatically deleted.
The BuildBot system does not take care of cleaning your experiment directory before or after each run. You may take precautions in your own run scripts for this.


h3. Signalling "SUCCESS" or "FAILURE"

Run scripts signal "SUCCESS" or "FAILURE" ("green" or "red") with the help of the script error status:
Ending the script with
<pre>
 exit 0
</pre>
(or omitting the explicit @exit@ call) means that the script has run successfully.
All other @exit@ values mean "FAILURE".

Two important remarks:

* The usual mechanism (on other BuildBot builders) requires that the scripts explicitly write the error status to a file @${PBS_O_WORKDIR}/${PBS_JOBNAME}.status@.
* The described @exit@ mechanism works only at DWD; for this a special handling has been implemented in @scripts/buildbot_scripts/runexp@: 
  Here, an alternative status file @${PBS_O_WORKDIR}/${PBS_JOBNAME}.status.2@ is automatically created.


h2. Admin Information


h3. User "buildbot"

A user "buildbot" is required. This user needs shell login, otherwise it is not possible to execute run scripts (though it is possible to submit them to the batch queue).
At DWD the user "buildbot" has the unusual home directory @/hpc/rhome/software/buildbot@. 
Make sure that the run scripts require a proper setting of the environment variable @$TMPDIR@.


h3. Access Rights

For security reasons, the user "buildbot" is not part of the group FE1 and therefore does not have read or write access to these parts of the file system.
Only the read access to the necessary input data is granted by setting access control lists (ACLs) for the subdirectory
<pre>
/lustre2/rwork0/routfor/test/icon/Checksuite_data
</pre>
by typing
<pre>
setfacl -d -R -m g:buildbot:rX /lustre2/rwork0/routfor/test/icon/Checksuite_data
setfacl  -R -m g:buildbot:rX /lustre2/rwork0/routfor/test/icon/Checksuite_data

getfacl /lustre2/rwork0/routfor/test/icon/Checksuite_data
</pre>


h3. Directory Structure for BuildBot at DWD

The BuildBot runs reside in the directory
<pre>
/lustre1/uwork/buildbot/icon/cray/default
</pre>
