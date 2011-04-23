#include <stdio.h>
#include <stdlib.h>

#include "config.h"

#ifdef HAVE_POSTGRESQL
#include "libpq-fe.h"
#endif

const int success = 0;
const int failed  = 1;

#if HAVE_POSTGRESQL == 1

static PGconn *connection;

int open_db_connection(char *connectionInfo)
{
  connection = PQconnectdb(connectionInfo);
  if (PQstatus(connection) != CONNECTION_OK)
    {
      fprintf(stderr, "Connection not established: %s\n", 
	      PQerrorMessage(connection));
      PQfinish(connection);
      return failed;
    }
  return success;
}

int begin_transaction_block(void)
{
  PGresult *res;

  res = PQexec(connection, "BEGIN");
  if (PQresultStatus(res) != PGRES_COMMAND_OK)
    {
      fprintf(stderr, "BEGIN command failed\n");
      PQclear(res);
      PQfinish(connection); 
      return failed;
    }
  PQclear(res);
  return success;
}

int commit_transaction_block(void)
{
  PGresult *res;

  res = PQexec(connection, "COMMIT");
  if (PQresultStatus(res) != PGRES_COMMAND_OK)
    {
      fprintf(stderr, "COMMIT command failed\n");
      PQclear(res);
      PQfinish(connection); 
      return failed;
    }
  PQclear(res);
  return success;
}

int insert_db_query(char *query)
{
  PGresult *res;

  res = PQexec(connection, query);
  if (PQresultStatus(res) != PGRES_COMMAND_OK)
    {
      fprintf(stderr, "Insertion of query failed: %s\n", 
	      PQerrorMessage(connection));
      PQclear(res);
      /* PQfinish(connection); */
      return failed;
    }
  PQclear(res);
  return success;
}

void close_db_connection(void)
{
  PQfinish(connection);
  return;
}

#else

int open_db_connection(void)
{
  return success;
}

int begin_transaction_block(void)
{
  return success;
}

int commit_transaction_block(void)
{
  return success;
}

int insert_db_query(char *experiment)
{
  return success;
}

void close_db_connection(void)
{
  return;
}

#endif
