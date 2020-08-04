#ifndef UTIL_ARITHMETIC_EXPR_H
#define UTIL_ARITHMETIC_EXPR_H


#define MAX_BUF_LEN     1024
#define NUM_FCT            9
#define MAX_NAME_LEN      32

typedef enum { VARIABLE, VALUE, OPERATOR, FUNCTION } expr_type;
typedef enum { EXP = 0, LOG = 1, SIN = 2, COS = 3, 
               MIN = 4, MAX = 5, IF  = 6, SQRT = 7,
               ERF = 8                              } fct_type;

extern const char * const fct_name[NUM_FCT];

/* Data type for postfix stack and queue (shunting yard algorithm) */
struct t_item {
  expr_type type;

  char      op;
  double    val;
  fct_type  fct;
  char      field[MAX_NAME_LEN];
};

/* Data type for stack or queue. */
struct t_list {
  int           size;
  struct t_item list[1024];
};



/* takes a string as input and generates a post-fix parse tree. */
int do_parse_infix(const char *in_parse_line, struct t_list* data);

/* accessor function for function name strings */
int get_fctname(const int id, char *string);

#endif /* UTIL_ARITHMETIC_EXPR_H */
