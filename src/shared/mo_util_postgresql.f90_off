MODULE mo_util_postgresql
  
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_DOUBLE, C_NULL_CHAR
  
  IMPLICIT NONE
  
  PRIVATE

  INTERFACE
    FUNCTION open_db_connection(connection_info) RESULT(iret) &
         BIND(C,NAME='open_db_connection')
      IMPORT :: C_CHAR,C_INT
      INTEGER(C_INT) :: iret
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: connection_info
    END FUNCTION open_db_connection
  END INTERFACE
  
  INTERFACE
    FUNCTION begin_transaction_block() RESULT(iret) &
         BIND(C,NAME='begin_transaction_block')
      IMPORT :: C_INT
      INTEGER(C_INT) :: iret
    END FUNCTION begin_transaction_block
  END INTERFACE

  INTERFACE
    FUNCTION commit_transaction_block() RESULT(iret) &
         BIND(C,NAME='commit_transaction_block')
      IMPORT :: C_INT
      INTEGER(C_INT) :: iret
    END FUNCTION commit_transaction_block
  END INTERFACE

  INTERFACE
    FUNCTION insert_db_query(query) RESULT(iret) &
         BIND(C,NAME='insert_dq_query')
      IMPORT :: C_CHAR,C_INT
      INTEGER(C_INT) :: iret
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in) :: query
    END FUNCTION insert_db_query
  END INTERFACE
  
  INTERFACE
    SUBROUTINE close_db_connection() BIND(C,NAME='close_db_connection')
    END SUBROUTINE close_db_connection
  END INTERFACE
  
  PUBLIC :: open_db_connection
  PUBLIC :: insert_query
  PUBLIC :: close_db_connection
  
  PUBLIC :: begin_transaction_block
  PUBLIC :: commit_transaction_block
  
CONTAINS
  
  FUNCTION insert_query(query) RESULT(iret)
    CHARACTER(len=*), INTENT(in) :: query 
    INTEGER(C_INT) :: iret
    iret = insert_db_query(TRIM(query)//C_NULL_CHAR)
  END FUNCTION insert_query
  
END MODULE mo_util_postgresql

