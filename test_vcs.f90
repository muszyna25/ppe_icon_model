PROGRAM test_vcs

  USE mo_util_vcs

  IMPLICIT NONE

  CHARACTER(len=256) :: repository, branch, revision
  INTEGER :: strlen

  strlen = LEN(repository)
  CALL util_repository_url(repository, strlen)
  PRINT *, repository(1:strlen)

  strlen = LEN(branch)
  CALL util_branch_name(branch, strlen)
  PRINT *, branch(1:strlen)

  strlen = LEN(branch)
  CALL util_revision_key(revision, strlen)
  PRINT *, revision(1:strlen)

END PROGRAM test_vcs
