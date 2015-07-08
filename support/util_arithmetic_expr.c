
#line 1 "util_arithmetic_expr.rl"
/* --------------------------------------------------------------------- *
 * Evaluation of basic arithmetic infix expressions based on a           *
 * Finite State Machine (FSM) and Dijkstra's shunting yard algorithm.    *
 * This file must be processed with Ragel to produce the final C code.   *
 *                                                                       *
 * 03/2015 : F. Prill, DWD                                               *
 * --------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util_arithmetic_expr.h"

const char*  const priorities = "(;><;+-;*/;^";
const char * const fct_name[NUM_FCT] = { "exp", "log", "sin", "cos", "min", "max", "if", "sqrt" };

/* All data needed to define the state of the Finite State Machine parser. */
struct t_parsedata {
  int  cs, act, buflen;                    /* FSM state, buffer length      */
  char *ts, *te;                           /* token start, end              */
  char buffer[MAX_BUF_LEN];                /* token string buffer           */
  struct t_list ostack, *queue;            /* operator stack, postfix queue */
};


/* Priority of an arithmetic operator */
int priority(char op) {
  const char* o = priorities;
  for (int val=0; o!='\0'; o++)
    if ((*o) == ';')     val++;
    else if ((*o) == op) return val;
  return -1;
}


/* --------------------------------------------------------------------- *
 * DEFINITION OF RAGEL FINITE STATE MACHINE "parse_infix"
 * --------------------------------------------------------------------- */


#line 158 "util_arithmetic_expr.rl"



#line 49 "util_arithmetic_expr.c"
static const char _parse_infix_actions[] = {
	0, 1, 0, 1, 1, 1, 2, 1, 
	3, 1, 4, 1, 5, 1, 6, 1, 
	7, 1, 8, 1, 9, 1, 10, 2, 
	1, 0
};

static const char _parse_infix_key_offsets[] = {
	0, 0, 24, 31, 39, 40, 41, 42, 
	43, 44, 45, 46, 48, 49, 50, 51, 
	53, 54, 55, 79, 82, 85, 88, 93, 
	99, 102, 107, 110, 113
};

static const char _parse_infix_trans_keys[] = {
	32, 44, 45, 47, 60, 62, 91, 94, 
	99, 101, 105, 108, 109, 112, 114, 115, 
	9, 13, 40, 41, 42, 43, 48, 57, 
	95, 48, 57, 65, 90, 97, 122, 93, 
	95, 48, 57, 65, 90, 97, 122, 111, 
	115, 120, 112, 102, 111, 103, 97, 105, 
	120, 110, 105, 105, 113, 114, 116, 32, 
	44, 45, 47, 60, 62, 91, 94, 99, 
	101, 105, 108, 109, 112, 114, 115, 9, 
	13, 40, 41, 42, 43, 48, 57, 32, 
	9, 13, 32, 9, 13, 32, 9, 13, 
	32, 9, 13, 48, 57, 32, 46, 9, 
	13, 48, 57, 32, 9, 13, 32, 9, 
	13, 48, 57, 32, 9, 13, 32, 9, 
	13, 32, 9, 13, 0
};

static const char _parse_infix_single_lengths[] = {
	0, 16, 1, 2, 1, 1, 1, 1, 
	1, 1, 1, 2, 1, 1, 1, 2, 
	1, 1, 16, 1, 1, 1, 1, 2, 
	1, 1, 1, 1, 1
};

static const char _parse_infix_range_lengths[] = {
	0, 4, 3, 3, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 4, 1, 1, 1, 2, 2, 
	1, 2, 1, 1, 1
};

static const unsigned char _parse_infix_index_offsets[] = {
	0, 0, 21, 26, 32, 34, 36, 38, 
	40, 42, 44, 46, 49, 51, 53, 55, 
	58, 60, 62, 83, 86, 89, 92, 96, 
	101, 104, 108, 111, 114
};

static const char _parse_infix_indicies[] = {
	0, 4, 5, 3, 3, 3, 7, 3, 
	8, 9, 10, 11, 12, 13, 14, 15, 
	0, 2, 3, 6, 1, 16, 16, 16, 
	16, 1, 17, 16, 16, 16, 16, 1, 
	18, 1, 19, 1, 20, 1, 19, 1, 
	19, 1, 21, 1, 19, 1, 22, 23, 
	1, 19, 1, 19, 1, 24, 1, 23, 
	25, 1, 26, 1, 19, 1, 0, 4, 
	5, 3, 3, 3, 7, 3, 8, 9, 
	10, 11, 12, 13, 14, 15, 0, 2, 
	3, 6, 1, 28, 28, 27, 30, 30, 
	29, 4, 4, 31, 30, 30, 32, 29, 
	34, 35, 34, 32, 33, 34, 34, 33, 
	34, 34, 35, 33, 17, 17, 36, 38, 
	38, 37, 40, 40, 39, 0
};

static const char _parse_infix_trans_targs[] = {
	1, 0, 19, 20, 21, 22, 23, 2, 
	4, 6, 8, 9, 11, 14, 28, 15, 
	3, 26, 5, 27, 7, 10, 12, 13, 
	28, 16, 17, 18, 19, 18, 20, 18, 
	23, 18, 24, 25, 18, 18, 27, 18, 
	28
};

static const char _parse_infix_trans_actions[] = {
	0, 0, 23, 23, 0, 23, 23, 3, 
	23, 23, 23, 23, 23, 23, 23, 23, 
	1, 0, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 15, 0, 13, 0, 17, 
	1, 9, 0, 1, 21, 11, 0, 19, 
	0
};

static const char _parse_infix_to_state_actions[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 5, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0
};

static const char _parse_infix_from_state_actions[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 7, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0
};

static const unsigned char _parse_infix_eof_trans[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 28, 30, 32, 30, 34, 
	34, 34, 37, 38, 40
};

static const int parse_infix_start = 18;
static const int parse_infix_first_final = 18;
static const int parse_infix_error = 0;

static const int parse_infix_en_main = 18;


#line 161 "util_arithmetic_expr.rl"


/* --------------------------------------------------------------------- */
int do_parse_infix(const char *in_parse_line, struct t_list *queue) {
  char *p0 = malloc(strlen (in_parse_line) + 1);     /* Space for length plus null */
  if (p0 == NULL) return -1;                         /* No memory                  */
  strcpy(p0, in_parse_line);                         /* Copy the characters        */

  char *p = p0, *pe = p0 + strlen(p0), *eof = pe ;   /* pointer to input state and end. */
  struct t_parsedata pdata;
  pdata.queue = queue;
  pdata.ostack.size = pdata.queue->size = 0;
  struct t_parsedata *data = &pdata;
  
  
#line 183 "util_arithmetic_expr.c"
	{
	 data->cs = parse_infix_start;
	 data->ts = 0;
	 data->te = 0;
	 data->act = 0;
	}

#line 177 "util_arithmetic_expr.rl"
  
#line 193 "util_arithmetic_expr.c"
	{
	int _klen;
	unsigned int _trans;
	const char *_acts;
	unsigned int _nacts;
	const char *_keys;

	if ( p == pe )
		goto _test_eof;
	if (  data->cs == 0 )
		goto _out;
_resume:
	_acts = _parse_infix_actions + _parse_infix_from_state_actions[ data->cs];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 ) {
		switch ( *_acts++ ) {
	case 3:
#line 1 "NONE"
	{ data->ts = p;}
	break;
#line 214 "util_arithmetic_expr.c"
		}
	}

	_keys = _parse_infix_trans_keys + _parse_infix_key_offsets[ data->cs];
	_trans = _parse_infix_index_offsets[ data->cs];

	_klen = _parse_infix_single_lengths[ data->cs];
	if ( _klen > 0 ) {
		const char *_lower = _keys;
		const char *_mid;
		const char *_upper = _keys + _klen - 1;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + ((_upper-_lower) >> 1);
			if ( (*p) < *_mid )
				_upper = _mid - 1;
			else if ( (*p) > *_mid )
				_lower = _mid + 1;
			else {
				_trans += (unsigned int)(_mid - _keys);
				goto _match;
			}
		}
		_keys += _klen;
		_trans += _klen;
	}

	_klen = _parse_infix_range_lengths[ data->cs];
	if ( _klen > 0 ) {
		const char *_lower = _keys;
		const char *_mid;
		const char *_upper = _keys + (_klen<<1) - 2;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + (((_upper-_lower) >> 1) & ~1);
			if ( (*p) < _mid[0] )
				_upper = _mid - 2;
			else if ( (*p) > _mid[1] )
				_lower = _mid + 2;
			else {
				_trans += (unsigned int)((_mid - _keys)>>1);
				goto _match;
			}
		}
		_trans += _klen;
	}

_match:
	_trans = _parse_infix_indicies[_trans];
_eof_trans:
	 data->cs = _parse_infix_trans_targs[_trans];

	if ( _parse_infix_trans_actions[_trans] == 0 )
		goto _again;

	_acts = _parse_infix_actions + _parse_infix_trans_actions[_trans];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 )
	{
		switch ( *_acts++ )
		{
	case 0:
#line 46 "util_arithmetic_expr.rl"
	{ if ((data->buflen < MAX_BUF_LEN) && ((*p)!='\n'))
   data->buffer[data->buflen++] = (*p);
 }
	break;
	case 1:
#line 50 "util_arithmetic_expr.rl"
	{ data->buflen = 0; }
	break;
	case 4:
#line 52 "util_arithmetic_expr.rl"
	{ data->te = p;p--;{ 
   data->queue->list[data->queue->size].type = VALUE;
   if (data->buflen < MAX_BUF_LEN) data->buffer[data->buflen] = 0;	
   sscanf(data->buffer, "%lf", &data->queue->list[data->queue->size++].val);
 }}
	break;
	case 5:
#line 68 "util_arithmetic_expr.rl"
	{ data->te = p;p--;{ 
   if (data->buflen < MAX_BUF_LEN) data->buffer[data->buflen] = 0;	
   const char* fct = data->buffer;
   struct t_item new_item =  { .type=FUNCTION };
   for (int i=0; i<NUM_FCT; i++)
     if (strcmp(fct, fct_name[i]) == 0)
       { new_item.fct = i; data->ostack.list[data->ostack.size++] = new_item; break; }
 }}
	break;
	case 6:
#line 86 "util_arithmetic_expr.rl"
	{ data->te = p;p--;{ 
  /*  While the stack is not empty and an operator is at the top and
      the operator at the top is higher priority that the item then
      pop the operator on the top of the stack, add the popped
      operator to the queue. */
  char this_op = data->buffer[0];
  while ((data->ostack.size > 0) && (data->ostack.list[data->ostack.size-1].type==OPERATOR) &&
	 (priority(this_op) < priority(data->ostack.list[data->ostack.size-1].op)))
    data->queue->list[data->queue->size++] = data->ostack.list[--data->ostack.size];
  struct t_item new_item =  { .type=OPERATOR, .op=this_op };
  data->ostack.list[data->ostack.size++] = new_item;
 }}
	break;
	case 7:
#line 99 "util_arithmetic_expr.rl"
	{ data->te = p;p--;{ 
  switch (data->buffer[0]) {
     case '(': {
         struct t_item new_item =  { .type=OPERATOR, .op='(' };
         data->ostack.list[data->ostack.size++] = new_item;
       }
       break;
  case ')':
       if (data->ostack.size > 0) {
         /* Until token at the top of the stack is a left parenthesis,
	    pop operators off the stack onto the output queue. */
         struct t_item top = data->ostack.list[--data->ostack.size];
         while ((data->ostack.size>0) && (top.type==OPERATOR) && (top.op != '(')) {
          data->queue->list[data->queue->size++] = top;
          top = data->ostack.list[--data->ostack.size];
         }
	 /* If the token at the top of the stack is a function token,
	    pop it onto the output queue. */
	 if (data->ostack.list[data->ostack.size-1].type==FUNCTION)
	   data->queue->list[data->queue->size++] = data->ostack.list[--data->ostack.size];
       } 
       break;        
  }
 }}
	break;
	case 8:
#line 77 "util_arithmetic_expr.rl"
	{ data->te = p;p--;{
  /* If the token is a function argument separator (i.e. a comma):
     Until the token at the top of the stack is a left parenthesis, pop
     operators off the stack onto the output queue. */
  while ((data->ostack.size>0) && (data->ostack.list[data->ostack.size-1].type==OPERATOR) && 
	 (data->ostack.list[data->ostack.size-1].op != '(')) 
    data->queue->list[data->queue->size++] = data->ostack.list[--data->ostack.size];
}}
	break;
	case 9:
#line 58 "util_arithmetic_expr.rl"
	{ data->te = p;p--;{ 
   if (data->buflen < MAX_BUF_LEN) data->buffer[data->buflen] = 0;	
   const char* constant = data->buffer;
   data->queue->list[data->queue->size].type = VALUE;
   if (strcmp(constant, "pi") == 0)
     data->queue->list[data->queue->size++].val = 4.*atan(1.);
   else if (strcmp(constant, "r") == 0)
     data->queue->list[data->queue->size++].val = 6.371229e6;
 }}
	break;
	case 10:
#line 124 "util_arithmetic_expr.rl"
	{ data->te = p;p--;{ 
   data->queue->list[data->queue->size].type = VARIABLE;
   if (data->buflen < MAX_BUF_LEN) data->buffer[data->buflen] = 0;	
   strcpy(data->queue->list[data->queue->size++].field, data->buffer);
 }}
	break;
#line 382 "util_arithmetic_expr.c"
		}
	}

_again:
	_acts = _parse_infix_actions + _parse_infix_to_state_actions[ data->cs];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 ) {
		switch ( *_acts++ ) {
	case 2:
#line 1 "NONE"
	{ data->ts = 0;}
	break;
#line 395 "util_arithmetic_expr.c"
		}
	}

	if (  data->cs == 0 )
		goto _out;
	if ( ++p != pe )
		goto _resume;
	_test_eof: {}
	if ( p == eof )
	{
	if ( _parse_infix_eof_trans[ data->cs] > 0 ) {
		_trans = _parse_infix_eof_trans[ data->cs] - 1;
		goto _eof_trans;
	}
	}

	_out: {}
	}

#line 179 "util_arithmetic_expr.rl"
  free(p0);
  if (data->cs == parse_infix_error) return 1;

  /* While the stack is not empty pop the stack, add to queue. */
  while (data->ostack.size > 0) 
    data->queue->list[data->queue->size++] = data->ostack.list[--data->ostack.size];
  return 0;
}


int get_fctname(const int id, char *string) {
  int ierr = 0;
  if ((id < 0) || (id >= NUM_FCT)) {
    ierr = 1;
  }
  else {
    strcpy(string, fct_name[id]);
  }
  return ierr;
}

