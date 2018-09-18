
#line 1 "usedep_parser.rl"
// --------------------------------------------------------------------------------
// Tokenizer for Fortran files; filters USE dependencies depending on 
// preprocessor #ifdef's, #define, and #undefine.
//
// Usage: usedep_parser  <root> "<defined symbols>" "<filenames>" [-v]
//
// Note:   This file must be processed with Ragel to produce the final C++ code.
//
// In more detail: This program contains two different parsers. First,
// a top-level scanner/parser separates #ifdef's and USE dependencies
// from the remaing content. A second-level scanner/parser then
// constructs an expression-syntax tree for each #ifdef condition.
//
// 11/2016: F. Prill, DWD
//
// --------------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>
#include <algorithm>
#include <memory> 
#include <map>
#include <iomanip>
#include <set> 
#include <sys/stat.h>


// --------------------------------------------------------------------------------
// Utility class defining a generic exception.
// --------------------------------------------------------------------------------

class ParserException : public std::exception
{
 private:
  std::string msg; // error message

 public:
  ParserException(const std::string msg) : msg(msg) {}
  virtual ~ParserException() throw() {}

  virtual const char* what() const throw()
  { return msg.c_str(); }
};


// --------------------------------------------------------------------------------
// Base class for ExpressionTokenizer and FileTokenizer
// --------------------------------------------------------------------------------

class Tokenizer
{
public:
  Tokenizer(const std::string& input) : buffer(input) 
  {
    // set up the buffer here
    p   = buffer.c_str();
    eof = pe = p + buffer.size();
  }

  std::string  buffer;                  //< string buffer
  const char  *ts, *pe;
protected:
  const char  *te, *p, *eof;
  int          act, cs, top, stack[1];  //< machine state
  std::string  cbuffer;                 //< buffer for current token
  
  const std::vector<std::string> PRIORITIES{"(", "||", "<", ">", "&&", "!"};

  // Priority of arithmetic operator, ie. pos in "PRIORITIES" array:
  int priority(std::string& op) {
    return (std::find(PRIORITIES.begin(), PRIORITIES.end(), op) - PRIORITIES.begin());
  }
};


// --------------------------------------------------------------------------------
// Expression Parser Definition: constructs an expression tree for IF-conditions
// --------------------------------------------------------------------------------


#line 132 "usedep_parser.rl"



#line 90 "usedep_parser.cpp"
static const char _ExpressionParser_actions[] = {
	0, 1, 1, 1, 2, 1, 3, 1, 
	6, 1, 7, 1, 8, 1, 9, 1, 
	10, 1, 11, 1, 12, 2, 0, 1, 
	2, 1, 0, 2, 4, 1, 3, 0, 
	1, 5
};

static const unsigned char _ExpressionParser_key_offsets[] = {
	0, 0, 2, 11, 19, 30, 34, 54, 
	57, 62, 64, 71, 80, 89, 98, 107, 
	116, 125, 136, 146, 149
};

static const char _ExpressionParser_trans_keys[] = {
	48, 57, 32, 40, 95, 9, 13, 65, 
	90, 97, 122, 32, 95, 9, 13, 65, 
	90, 97, 122, 32, 41, 95, 9, 13, 
	48, 57, 65, 90, 97, 122, 32, 41, 
	9, 13, 32, 33, 38, 45, 60, 62, 
	68, 95, 100, 124, 9, 13, 40, 41, 
	48, 57, 65, 90, 97, 122, 32, 9, 
	13, 33, 38, 60, 62, 124, 48, 57, 
	95, 48, 57, 65, 90, 97, 122, 69, 
	95, 101, 48, 57, 65, 90, 97, 122, 
	70, 95, 102, 48, 57, 65, 90, 97, 
	122, 73, 95, 105, 48, 57, 65, 90, 
	97, 122, 78, 95, 110, 48, 57, 65, 
	90, 97, 122, 69, 95, 101, 48, 57, 
	65, 90, 97, 122, 68, 95, 100, 48, 
	57, 65, 90, 97, 122, 32, 40, 95, 
	9, 13, 48, 57, 65, 90, 97, 122, 
	32, 95, 9, 13, 48, 57, 65, 90, 
	97, 122, 32, 9, 13, 32, 95, 9, 
	13, 48, 57, 65, 90, 97, 122, 0
};

static const char _ExpressionParser_single_lengths[] = {
	0, 0, 3, 2, 3, 2, 10, 1, 
	5, 0, 1, 3, 3, 3, 3, 3, 
	3, 3, 2, 1, 2
};

static const char _ExpressionParser_range_lengths[] = {
	0, 1, 3, 3, 4, 1, 5, 1, 
	0, 1, 3, 3, 3, 3, 3, 3, 
	3, 4, 4, 1, 4
};

static const unsigned char _ExpressionParser_index_offsets[] = {
	0, 0, 2, 9, 15, 23, 27, 43, 
	46, 52, 54, 59, 66, 73, 80, 87, 
	94, 101, 109, 116, 119
};

static const char _ExpressionParser_indicies[] = {
	0, 1, 3, 4, 5, 3, 5, 5, 
	2, 4, 6, 4, 6, 6, 2, 7, 
	8, 9, 7, 9, 9, 9, 2, 7, 
	8, 7, 2, 10, 11, 11, 13, 11, 
	11, 16, 15, 16, 11, 10, 12, 14, 
	15, 15, 1, 10, 10, 17, 19, 19, 
	19, 19, 19, 18, 0, 20, 22, 22, 
	22, 22, 21, 23, 22, 23, 22, 22, 
	22, 21, 24, 22, 24, 22, 22, 22, 
	21, 25, 22, 25, 22, 22, 22, 21, 
	26, 22, 26, 22, 22, 22, 21, 27, 
	22, 27, 22, 22, 22, 21, 28, 22, 
	28, 22, 22, 22, 21, 3, 4, 29, 
	3, 22, 29, 29, 21, 31, 32, 31, 
	32, 32, 32, 30, 31, 31, 30, 31, 
	33, 31, 33, 33, 33, 21, 0
};

static const char _ExpressionParser_trans_targs[] = {
	9, 0, 6, 2, 3, 18, 4, 5, 
	6, 4, 7, 8, 6, 1, 9, 10, 
	11, 6, 6, 8, 6, 6, 10, 12, 
	13, 14, 15, 16, 17, 20, 6, 19, 
	18, 20
};

static const char _ExpressionParser_trans_actions[] = {
	1, 0, 19, 0, 0, 21, 21, 0, 
	7, 1, 0, 21, 30, 21, 21, 21, 
	21, 17, 15, 1, 9, 11, 1, 1, 
	1, 1, 1, 1, 27, 24, 13, 0, 
	1, 1
};

static const char _ExpressionParser_to_state_actions[] = {
	0, 0, 0, 0, 0, 0, 3, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0
};

static const char _ExpressionParser_from_state_actions[] = {
	0, 0, 0, 0, 0, 0, 5, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0
};

static const unsigned char _ExpressionParser_eof_trans[] = {
	0, 0, 3, 3, 3, 3, 0, 18, 
	19, 21, 22, 22, 22, 22, 22, 22, 
	22, 22, 31, 31, 22
};

static const int ExpressionParser_start = 6;
static const int ExpressionParser_error = 0;

static const int ExpressionParser_en_main = 6;


#line 135 "usedep_parser.rl"


class ExpressionTokenizer : public Tokenizer
{
public:
  typedef enum { NONE, OPERATOR, VALUE, VARIABLE } expr_type;
  
  /* Data type for postfix stack and queue (shunting yard algorithm) */
  class t_item {
  public:
    t_item(expr_type type = NONE, std::string val = "")
      : type(type), val(val) {}
    expr_type    type;
    std::string  val;
  };

  ExpressionTokenizer(const std::string& input) : Tokenizer(input) 
  {
    
#line 225 "usedep_parser.cpp"
	{
	cs = ExpressionParser_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 154 "usedep_parser.rl"
    
#line 235 "usedep_parser.cpp"
	{
	int _klen;
	unsigned int _trans;
	const char *_acts;
	unsigned int _nacts;
	const char *_keys;

	if ( p == pe )
		goto _test_eof;
	if ( cs == 0 )
		goto _out;
_resume:
	_acts = _ExpressionParser_actions + _ExpressionParser_from_state_actions[cs];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 ) {
		switch ( *_acts++ ) {
	case 3:
#line 1 "NONE"
	{ts = p;}
	break;
#line 256 "usedep_parser.cpp"
		}
	}

	_keys = _ExpressionParser_trans_keys + _ExpressionParser_key_offsets[cs];
	_trans = _ExpressionParser_index_offsets[cs];

	_klen = _ExpressionParser_single_lengths[cs];
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

	_klen = _ExpressionParser_range_lengths[cs];
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
	_trans = _ExpressionParser_indicies[_trans];
_eof_trans:
	cs = _ExpressionParser_trans_targs[_trans];

	if ( _ExpressionParser_trans_actions[_trans] == 0 )
		goto _again;

	_acts = _ExpressionParser_actions + _ExpressionParser_trans_actions[_trans];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 )
	{
		switch ( *_acts++ )
		{
	case 0:
#line 85 "usedep_parser.rl"
	{ cbuffer = "";  }
	break;
	case 1:
#line 86 "usedep_parser.rl"
	{ if ((*p)!='\n') cbuffer += (*p); }
	break;
	case 4:
#line 1 "NONE"
	{te = p+1;}
	break;
	case 5:
#line 101 "usedep_parser.rl"
	{te = p+1;{
      if (cbuffer=="(") 
	ostack.push_back( t_item(OPERATOR, "(") );
      else if ((cbuffer==")") && (!ostack.empty())) {
	/* Until token at the top of the stack is a left parenthesis,
	   pop operators off the stack onto the output queue. */
	t_item top = ostack.back(); ostack.pop_back();
	while ((top.type==OPERATOR) && (top.val != "(")) {
	  queue.push_back(top);
          if (ostack.empty())  throw(ParserException("No matching opening parenthesis!"));
	  top = ostack.back(); ostack.pop_back();
	}
      }
    }}
	break;
	case 6:
#line 126 "usedep_parser.rl"
	{te = p+1;{ queue.push_back( t_item(VARIABLE, cbuffer) ); }}
	break;
	case 7:
#line 124 "usedep_parser.rl"
	{te = p;p--;{ queue.push_back( t_item(VALUE,    cbuffer) ); }}
	break;
	case 8:
#line 125 "usedep_parser.rl"
	{te = p;p--;{ queue.push_back( t_item(VARIABLE, cbuffer) ); }}
	break;
	case 9:
#line 127 "usedep_parser.rl"
	{te = p;p--;{ queue.push_back( t_item(VARIABLE, cbuffer) ); }}
	break;
	case 10:
#line 89 "usedep_parser.rl"
	{te = p;p--;{ 
      /*  While the stack is not empty and an operator is at the top
          and the operator at the top is higher priority than the item
          then pop the operator on the top of the stack, add the
          popped operator to the queue. */
      int pr = priority(cbuffer);
      while ((!ostack.empty()) && (ostack.back().type==OPERATOR) &&
	     (pr < priority(ostack.back().val)))
        { queue.push_back(ostack.back()); ostack.pop_back(); }
      ostack.push_back( t_item(OPERATOR, cbuffer) );
     }}
	break;
	case 11:
#line 129 "usedep_parser.rl"
	{te = p;p--;}
	break;
	case 12:
#line 125 "usedep_parser.rl"
	{{p = ((te))-1;}{ queue.push_back( t_item(VARIABLE, cbuffer) ); }}
	break;
#line 389 "usedep_parser.cpp"
		}
	}

_again:
	_acts = _ExpressionParser_actions + _ExpressionParser_to_state_actions[cs];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 ) {
		switch ( *_acts++ ) {
	case 2:
#line 1 "NONE"
	{ts = 0;}
	break;
#line 402 "usedep_parser.cpp"
		}
	}

	if ( cs == 0 )
		goto _out;
	if ( ++p != pe )
		goto _resume;
	_test_eof: {}
	if ( p == eof )
	{
	if ( _ExpressionParser_eof_trans[cs] > 0 ) {
		_trans = _ExpressionParser_eof_trans[cs] - 1;
		goto _eof_trans;
	}
	}

	_out: {}
	}

#line 155 "usedep_parser.rl"

    if (cs == ExpressionParser_error) 
      throw(ParserException("Parsing of condition failed for \""+buffer+"\""));
    
    // while the stack is not empty pop the stack, add to queue.
    while (!ostack.empty()) { 
      if (ostack.back().val == "(")   throw(ParserException("Unbalanced parentheses!"));
      queue.push_back(ostack.back()); ostack.pop_back(); 
    }
  }

  void negate() 
  { queue.push_back( t_item(OPERATOR, "!") ); } // negate the whole condition

  void print_stack() {
    std::cout << "Queue:" << std::endl;
    for (unsigned int i=0; i<queue.size(); ++i)
      std::cout << "i=" << i << ": type=" << queue[i].type << " , " << queue[i].val << std::endl;
  }

  int evaluate(std::map<std::string,int> defined_symbols) 
  {
    std::vector<int> rstack;
    for (unsigned int i=0; i<queue.size(); i++)
      {
        switch (queue[i].type) {
        case VALUE:
          try {
            rstack.push_back(std::stoi(queue[i].val));
          } 
          catch (std::exception& e) {
            throw(ParserException("Number conversion failed for \""+buffer+"\")"));
          }
          break;
        case VARIABLE:
          {
            if (defined_symbols.find(queue[i].val) != defined_symbols.end())
              rstack.push_back(defined_symbols[queue[i].val]);
            else
              rstack.push_back(0); // (undefined symbol)
          }
          break;
        case OPERATOR:
          {
            double op2 = rstack.back(); rstack.pop_back(); 
            if (queue[i].val ==  "!")
              { rstack.push_back( !op2 ); } // unary operator
            else if (!rstack.empty()) {
              double op1 = rstack.back(); rstack.pop_back(); 
              if      (queue[i].val == "&&")
                { rstack.push_back( (op1 > 0) && (op2 > 0) ); }
              else if (queue[i].val == "||")
                { rstack.push_back( (op1 > 0) || (op2 > 0) ); }
              else if (queue[i].val ==  ">")
                { rstack.push_back( op1 > op2 ); }
              else if (queue[i].val ==  "<")
                { rstack.push_back( op1 < op2 ); }
              else
                throw(ParserException("Internal error: unknown operator (evaluation for \""+buffer+"\")"));
            }
            else
              throw(ParserException("Internal error: not enough operands (evaluation for \""+buffer+"\")"));
          }
          break;
        case NONE:
          throw(ParserException("Internal error (evaluation for \""+buffer+"\")"));
          break;
        }
      }
    if (rstack.size() != 1)  throw(ParserException("Evaluation of stack failed for \""+buffer+"\""));
    return rstack[0];
  }

private:
  std::vector<t_item> ostack, queue;  //< operator stack, postfix queue
};


// --------------------------------------------------------------------------------
// FileTokenizer: separates #ifdef's and USE dependencies from the remaing content
// --------------------------------------------------------------------------------


#line 297 "usedep_parser.rl"



#line 510 "usedep_parser.cpp"
static const char _TokenParser_actions[] = {
	0, 1, 1, 1, 2, 1, 5, 1, 
	13, 1, 14, 1, 15, 1, 16, 1, 
	17, 1, 18, 1, 19, 1, 20, 1, 
	24, 1, 25, 1, 26, 1, 27, 1, 
	28, 1, 29, 1, 30, 1, 31, 2, 
	0, 21, 2, 0, 22, 2, 0, 23, 
	2, 1, 2, 2, 3, 4, 2, 6, 
	2, 2, 6, 7, 2, 6, 11, 2, 
	6, 12, 3, 6, 1, 2, 3, 6, 
	7, 1, 4, 6, 1, 2, 8, 4, 
	6, 1, 2, 9, 4, 6, 1, 2, 
	10
};

static const char _TokenParser_cond_offsets[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 1, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 3, 4, 5, 6, 
	6, 6, 6, 6, 6, 6, 6, 6, 
	6, 6, 6, 6, 6, 6, 6, 6, 
	6, 6, 6, 6
};

static const char _TokenParser_cond_lengths[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 1, 1, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 1, 1, 1, 1, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0
};

static const short _TokenParser_cond_keys[] = {
	101, 101, 101, 101, 110, 110, 100, 100, 
	105, 105, 102, 102, 0
};

static const char _TokenParser_cond_spaces[] = {
	0, 0, 0, 0, 0, 0, 0
};

static const short _TokenParser_key_offsets[] = {
	0, 0, 5, 6, 7, 8, 10, 12, 
	16, 18, 20, 22, 24, 26, 29, 37, 
	46, 48, 50, 52, 54, 56, 58, 61, 
	69, 78, 80, 82, 84, 87, 95, 108, 
	114, 115, 119, 126, 133, 142, 149, 156, 
	163, 170, 177, 185, 198, 211, 217, 224, 
	231, 238, 245, 252, 260, 273, 286, 292, 
	299, 306, 314, 327, 343, 352, 353, 361, 
	366, 368, 370, 376, 382, 383, 384, 385, 
	386, 387, 390, 398, 404, 411, 414, 418, 
	420, 422, 423, 424, 442, 443, 444, 445, 
	463, 464, 465, 466, 467, 468, 473, 477, 
	487, 488, 489, 490, 491, 494, 502, 503, 
	505, 506, 524, 525, 527, 528, 529, 530, 
	532, 537, 551, 565, 574, 580, 584, 592, 
	600, 608, 609, 610, 618, 619, 637, 653, 
	671, 687, 694, 712
};

static const short _TokenParser_trans_keys[] = {
	10, 33, 34, 35, 39, 10, 34, 39, 
	78, 110, 68, 100, 9, 32, 77, 109, 
	79, 111, 68, 100, 85, 117, 76, 108, 
	69, 101, 32, 9, 13, 32, 95, 9, 
	13, 65, 90, 97, 122, 10, 32, 95, 
	48, 57, 65, 90, 97, 122, 10, 32, 
	79, 111, 68, 100, 85, 117, 76, 108, 
	69, 101, 32, 9, 13, 32, 95, 9, 
	13, 65, 90, 97, 122, 10, 32, 95, 
	48, 57, 65, 90, 97, 122, 10, 32, 
	83, 115, 69, 101, 32, 9, 13, 32, 
	95, 9, 13, 65, 90, 97, 122, 10, 
	32, 33, 44, 95, 9, 13, 48, 57, 
	65, 90, 97, 122, 10, 32, 33, 44, 
	9, 13, 10, 10, 32, 9, 13, 10, 
	33, 34, 35, 39, 78, 110, 10, 33, 
	34, 35, 39, 68, 100, 9, 10, 32, 
	33, 34, 35, 39, 77, 109, 10, 33, 
	34, 35, 39, 79, 111, 10, 33, 34, 
	35, 39, 68, 100, 10, 33, 34, 35, 
	39, 85, 117, 10, 33, 34, 35, 39, 
	76, 108, 10, 33, 34, 35, 39, 69, 
	101, 10, 32, 33, 34, 35, 39, 9, 
	13, 10, 32, 33, 34, 35, 39, 95, 
	9, 13, 65, 90, 97, 122, 10, 32, 
	33, 34, 35, 39, 95, 48, 57, 65, 
	90, 97, 122, 10, 32, 33, 34, 35, 
	39, 10, 33, 34, 35, 39, 79, 111, 
	10, 33, 34, 35, 39, 68, 100, 10, 
	33, 34, 35, 39, 85, 117, 10, 33, 
	34, 35, 39, 76, 108, 10, 33, 34, 
	35, 39, 69, 101, 10, 32, 33, 34, 
	35, 39, 9, 13, 10, 32, 33, 34, 
	35, 39, 95, 9, 13, 65, 90, 97, 
	122, 10, 32, 33, 34, 35, 39, 95, 
	48, 57, 65, 90, 97, 122, 10, 32, 
	33, 34, 35, 39, 10, 33, 34, 35, 
	39, 83, 115, 10, 33, 34, 35, 39, 
	69, 101, 10, 32, 33, 34, 35, 39, 
	9, 13, 10, 32, 33, 34, 35, 39, 
	95, 9, 13, 65, 90, 97, 122, 10, 
	32, 33, 34, 35, 39, 44, 95, 9, 
	13, 48, 57, 65, 90, 97, 122, 10, 
	32, 33, 34, 35, 39, 44, 9, 13, 
	10, 10, 32, 33, 34, 35, 39, 9, 
	13, 10, 33, 34, 35, 39, 10, 34, 
	10, 39, 32, 100, 105, 117, 357, 613, 
	32, 100, 105, 117, 357, 613, 101, 102, 
	105, 110, 101, 32, 9, 13, 32, 95, 
	9, 13, 65, 90, 97, 122, 48, 57, 
	65, 90, 97, 122, 41, 48, 57, 65, 
	90, 97, 122, 32, 9, 13, 32, 35, 
	9, 13, 102, 110, 100, 110, 101, 102, 
	32, 33, 38, 45, 60, 62, 95, 124, 
	9, 13, 40, 41, 48, 57, 65, 90, 
	97, 122, 100, 101, 102, 32, 33, 38, 
	45, 60, 62, 95, 124, 9, 13, 40, 
	41, 48, 57, 65, 90, 97, 122, 99, 
	108, 117, 100, 101, 32, 34, 39, 9, 
	13, 65, 90, 97, 122, 34, 39, 46, 
	95, 48, 57, 65, 90, 97, 122, 110, 
	100, 101, 102, 32, 9, 13, 32, 95, 
	9, 13, 65, 90, 97, 122, 108, 105, 
	115, 102, 32, 33, 38, 45, 60, 62, 
	95, 124, 9, 13, 40, 41, 48, 57, 
	65, 90, 97, 122, 101, 108, 622, 612, 
	617, 614, 102, 110, 32, 100, 110, 9, 
	13, 10, 32, 33, 34, 35, 39, 69, 
	77, 85, 101, 109, 117, 9, 13, 10, 
	32, 33, 34, 35, 39, 69, 77, 85, 
	101, 109, 117, 9, 13, 32, 69, 77, 
	85, 101, 109, 117, 9, 13, 10, 32, 
	33, 44, 9, 13, 10, 32, 9, 13, 
	32, 95, 9, 13, 65, 90, 97, 122, 
	32, 95, 9, 13, 65, 90, 97, 122, 
	32, 95, 9, 13, 65, 90, 97, 122, 
	34, 39, 40, 95, 48, 57, 65, 90, 
	97, 122, 32, 32, 33, 38, 45, 60, 
	62, 95, 124, 9, 13, 40, 41, 48, 
	57, 65, 90, 97, 122, 38, 45, 60, 
	62, 95, 124, 32, 33, 40, 41, 48, 
	57, 65, 90, 97, 122, 32, 33, 38, 
	45, 60, 62, 95, 124, 9, 13, 40, 
	41, 48, 57, 65, 90, 97, 122, 38, 
	45, 60, 62, 95, 124, 32, 33, 40, 
	41, 48, 57, 65, 90, 97, 122, 95, 
	48, 57, 65, 90, 97, 122, 32, 33, 
	38, 45, 60, 62, 95, 124, 9, 13, 
	40, 41, 48, 57, 65, 90, 97, 122, 
	38, 45, 60, 62, 95, 124, 32, 33, 
	40, 41, 48, 57, 65, 90, 97, 122, 
	0
};

static const char _TokenParser_single_lengths[] = {
	0, 5, 1, 1, 1, 2, 2, 4, 
	2, 2, 2, 2, 2, 1, 2, 3, 
	2, 2, 2, 2, 2, 2, 1, 2, 
	3, 2, 2, 2, 1, 2, 5, 4, 
	1, 2, 7, 7, 9, 7, 7, 7, 
	7, 7, 6, 7, 7, 6, 7, 7, 
	7, 7, 7, 6, 7, 7, 6, 7, 
	7, 6, 7, 8, 7, 1, 6, 5, 
	2, 2, 6, 6, 1, 1, 1, 1, 
	1, 1, 2, 0, 1, 1, 2, 2, 
	2, 1, 1, 8, 1, 1, 1, 8, 
	1, 1, 1, 1, 1, 3, 0, 4, 
	1, 1, 1, 1, 1, 2, 1, 2, 
	1, 8, 1, 2, 1, 1, 1, 2, 
	3, 12, 12, 7, 4, 2, 2, 2, 
	2, 1, 1, 2, 1, 8, 6, 8, 
	6, 1, 8, 6
};

static const char _TokenParser_range_lengths[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 1, 3, 3, 
	0, 0, 0, 0, 0, 0, 1, 3, 
	3, 0, 0, 0, 1, 3, 4, 1, 
	0, 1, 0, 0, 0, 0, 0, 0, 
	0, 0, 1, 3, 3, 0, 0, 0, 
	0, 0, 0, 1, 3, 3, 0, 0, 
	0, 1, 3, 4, 1, 0, 1, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 1, 3, 3, 3, 1, 1, 0, 
	0, 0, 0, 5, 0, 0, 0, 5, 
	0, 0, 0, 0, 0, 1, 2, 3, 
	0, 0, 0, 0, 1, 3, 0, 0, 
	0, 5, 0, 0, 0, 0, 0, 0, 
	1, 1, 1, 1, 1, 1, 3, 3, 
	3, 0, 0, 3, 0, 5, 5, 5, 
	5, 3, 5, 5
};

static const short _TokenParser_index_offsets[] = {
	0, 0, 6, 8, 10, 12, 15, 18, 
	23, 26, 29, 32, 35, 38, 41, 47, 
	54, 57, 60, 63, 66, 69, 72, 75, 
	81, 88, 91, 94, 97, 100, 106, 116, 
	122, 124, 128, 136, 144, 154, 162, 170, 
	178, 186, 194, 202, 213, 224, 231, 239, 
	247, 255, 263, 271, 279, 290, 301, 308, 
	316, 324, 332, 343, 356, 365, 367, 375, 
	381, 384, 387, 394, 401, 403, 405, 407, 
	409, 411, 414, 420, 424, 429, 432, 436, 
	439, 442, 444, 446, 460, 462, 464, 466, 
	480, 482, 484, 486, 488, 490, 495, 498, 
	506, 508, 510, 512, 514, 517, 523, 525, 
	528, 530, 544, 546, 549, 551, 553, 555, 
	558, 563, 577, 591, 600, 606, 610, 616, 
	622, 628, 630, 632, 638, 640, 654, 666, 
	680, 692, 697, 711
};

static const unsigned char _TokenParser_indicies[] = {
	2, 3, 4, 0, 5, 1, 2, 3, 
	1, 4, 1, 5, 7, 7, 6, 8, 
	8, 6, 8, 8, 9, 9, 6, 10, 
	10, 6, 11, 11, 6, 12, 12, 6, 
	13, 13, 6, 14, 14, 6, 15, 15, 
	6, 15, 16, 15, 16, 16, 0, 17, 
	18, 19, 19, 19, 19, 0, 17, 18, 
	0, 20, 20, 6, 21, 21, 6, 22, 
	22, 6, 23, 23, 6, 24, 24, 6, 
	25, 25, 6, 25, 26, 25, 26, 26, 
	0, 27, 28, 29, 29, 29, 29, 0, 
	27, 28, 0, 30, 30, 6, 31, 31, 
	6, 32, 32, 6, 32, 33, 32, 33, 
	33, 0, 35, 34, 36, 37, 38, 34, 
	38, 38, 38, 0, 35, 34, 36, 37, 
	34, 0, 39, 36, 40, 37, 37, 36, 
	2, 3, 4, 0, 5, 41, 41, 1, 
	2, 3, 4, 0, 5, 42, 42, 1, 
	42, 2, 42, 3, 4, 0, 5, 43, 
	43, 1, 2, 3, 4, 0, 5, 44, 
	44, 1, 2, 3, 4, 0, 5, 45, 
	45, 1, 2, 3, 4, 0, 5, 46, 
	46, 1, 2, 3, 4, 0, 5, 47, 
	47, 1, 2, 3, 4, 0, 5, 48, 
	48, 1, 50, 49, 3, 4, 0, 5, 
	49, 1, 50, 49, 3, 4, 0, 5, 
	51, 49, 51, 51, 1, 17, 52, 3, 
	4, 0, 5, 53, 53, 53, 53, 1, 
	17, 52, 3, 4, 0, 5, 1, 2, 
	3, 4, 0, 5, 54, 54, 1, 2, 
	3, 4, 0, 5, 55, 55, 1, 2, 
	3, 4, 0, 5, 56, 56, 1, 2, 
	3, 4, 0, 5, 57, 57, 1, 2, 
	3, 4, 0, 5, 58, 58, 1, 60, 
	59, 3, 4, 0, 5, 59, 1, 60, 
	59, 3, 4, 0, 5, 61, 59, 61, 
	61, 1, 27, 62, 3, 4, 0, 5, 
	63, 63, 63, 63, 1, 27, 62, 3, 
	4, 0, 5, 1, 2, 3, 4, 0, 
	5, 64, 64, 1, 2, 3, 4, 0, 
	5, 65, 65, 1, 67, 66, 3, 4, 
	0, 5, 66, 1, 67, 66, 3, 4, 
	0, 5, 68, 66, 68, 68, 1, 35, 
	69, 70, 4, 0, 5, 71, 72, 69, 
	72, 72, 72, 1, 35, 69, 70, 4, 
	0, 5, 71, 69, 1, 39, 70, 40, 
	71, 70, 74, 36, 75, 71, 73, 39, 
	70, 74, 36, 75, 73, 76, 73, 74, 
	77, 73, 75, 78, 80, 81, 82, 83, 
	84, 79, 78, 80, 85, 82, 83, 84, 
	79, 86, 79, 87, 79, 88, 79, 89, 
	79, 90, 79, 91, 91, 79, 91, 92, 
	91, 92, 92, 79, 94, 94, 94, 93, 
	95, 94, 94, 94, 93, 96, 96, 93, 
	96, 97, 96, 93, 98, 99, 79, 100, 
	101, 79, 102, 79, 103, 79, 104, 105, 
	105, 105, 105, 105, 105, 105, 103, 105, 
	105, 105, 105, 0, 106, 79, 107, 79, 
	108, 79, 109, 110, 110, 110, 110, 110, 
	110, 110, 108, 110, 110, 110, 110, 0, 
	111, 79, 112, 79, 113, 79, 114, 79, 
	115, 79, 115, 116, 116, 115, 79, 117, 
	117, 79, 118, 118, 119, 119, 119, 119, 
	119, 79, 120, 79, 121, 79, 122, 79, 
	123, 79, 124, 124, 79, 124, 125, 124, 
	125, 125, 79, 126, 79, 127, 128, 79, 
	129, 79, 130, 131, 131, 131, 131, 131, 
	131, 131, 129, 131, 131, 131, 131, 0, 
	132, 79, 126, 133, 79, 134, 79, 135, 
	79, 136, 79, 137, 99, 79, 103, 100, 
	101, 103, 79, 140, 139, 141, 142, 143, 
	144, 145, 146, 147, 145, 146, 147, 139, 
	138, 150, 149, 3, 4, 148, 5, 151, 
	152, 153, 151, 152, 153, 149, 1, 150, 
	154, 155, 156, 154, 155, 156, 150, 148, 
	35, 34, 36, 37, 34, 157, 40, 37, 
	37, 36, 15, 16, 15, 16, 16, 158, 
	25, 26, 25, 26, 26, 158, 32, 33, 
	32, 33, 33, 158, 1, 4, 1, 5, 
	160, 161, 161, 161, 161, 159, 97, 162, 
	104, 105, 105, 105, 105, 105, 105, 105, 
	103, 105, 105, 105, 105, 163, 164, 164, 
	164, 164, 164, 164, 164, 164, 164, 164, 
	164, 163, 109, 110, 110, 110, 110, 110, 
	110, 110, 108, 110, 110, 110, 110, 165, 
	166, 166, 166, 166, 166, 166, 166, 166, 
	166, 166, 166, 165, 168, 168, 168, 168, 
	167, 130, 131, 131, 131, 131, 131, 131, 
	131, 129, 131, 131, 131, 131, 169, 170, 
	170, 170, 170, 170, 170, 170, 170, 170, 
	170, 170, 169, 0
};

static const unsigned char _TokenParser_trans_targs[] = {
	113, 1, 113, 2, 3, 4, 113, 6, 
	7, 8, 9, 10, 11, 12, 13, 14, 
	15, 113, 16, 15, 18, 19, 20, 21, 
	22, 23, 24, 113, 25, 24, 27, 28, 
	29, 30, 31, 116, 32, 33, 30, 113, 
	117, 35, 36, 37, 38, 39, 40, 41, 
	42, 43, 118, 44, 45, 44, 47, 48, 
	49, 50, 51, 52, 119, 53, 54, 53, 
	56, 57, 58, 120, 59, 60, 61, 62, 
	59, 63, 64, 65, 121, 122, 67, 0, 
	68, 111, 96, 102, 107, 79, 69, 70, 
	71, 72, 73, 74, 123, 113, 76, 77, 
	78, 124, 80, 88, 81, 84, 82, 83, 
	125, 126, 85, 86, 87, 127, 128, 89, 
	90, 91, 92, 93, 94, 95, 113, 95, 
	97, 98, 99, 100, 101, 129, 103, 104, 
	106, 105, 130, 131, 113, 108, 109, 110, 
	113, 112, 1, 114, 115, 2, 3, 66, 
	4, 34, 46, 55, 113, 114, 115, 34, 
	46, 55, 5, 17, 26, 113, 113, 113, 
	75, 123, 113, 113, 126, 113, 128, 113, 
	129, 113, 131
};

static const char _TokenParser_trans_actions[] = {
	37, 0, 19, 0, 0, 0, 33, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	48, 17, 0, 3, 0, 0, 0, 0, 
	0, 0, 48, 15, 0, 3, 0, 0, 
	0, 48, 0, 60, 0, 0, 3, 11, 
	60, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 63, 48, 0, 3, 0, 0, 
	0, 0, 0, 0, 63, 48, 0, 3, 
	0, 0, 0, 63, 48, 0, 0, 0, 
	3, 0, 0, 0, 60, 60, 0, 0, 
	0, 0, 0, 1, 1, 0, 0, 0, 
	0, 0, 0, 0, 66, 35, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	74, 48, 0, 0, 0, 79, 48, 0, 
	0, 0, 0, 0, 0, 48, 13, 3, 
	0, 0, 0, 0, 0, 48, 0, 0, 
	0, 0, 84, 48, 7, 0, 0, 0, 
	9, 0, 1, 70, 70, 1, 1, 1, 
	1, 1, 1, 1, 21, 57, 57, 0, 
	0, 0, 0, 0, 0, 23, 29, 25, 
	0, 54, 31, 39, 3, 42, 3, 27, 
	3, 45, 3
};

static const char _TokenParser_to_state_actions[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 51, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0
};

static const char _TokenParser_from_state_actions[] = {
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 5, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0
};

static const short _TokenParser_eof_trans[] = {
	0, 1, 1, 1, 1, 7, 7, 7, 
	7, 7, 7, 7, 7, 7, 1, 1, 
	1, 7, 7, 7, 7, 7, 7, 1, 
	1, 1, 7, 7, 7, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 
	1, 1, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 94, 94, 94, 94, 0, 
	0, 0, 0, 1, 0, 0, 0, 1, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 
	0, 1, 0, 0, 0, 0, 0, 0, 
	0, 0, 149, 149, 158, 158, 159, 159, 
	159, 158, 158, 160, 163, 164, 164, 166, 
	166, 168, 170, 170
};

static const int TokenParser_start = 113;
static const int TokenParser_first_final = 113;
static const int TokenParser_error = 0;

static const int TokenParser_en_main = 113;


#line 300 "usedep_parser.rl"


class FileTokenizer : public Tokenizer
{
public:
    
    struct Token
    { 
        typedef std::shared_ptr<Token> ptr;

        enum  TokenType                         
          { End,   IfThen,   Else,   ElseIf,   EndIf,   Usedep,   Code,   Incldep,   
            Module,   EndModule,   Define,   Undefine,   None };
        const std::vector<std::string> TokenName
          {"End", "IfThen", "Else", "ElseIf", "EndIf", "Usedep", "Code", "Incldep", 
           "Module", "EndModule", "Define", "Undefine", "None"};
    
        TokenType                             type;
        std::string                           value;
        std::unique_ptr<ExpressionTokenizer>  condition;
    };

    class TreeNode
    {
      public:
        typedef std::shared_ptr<TreeNode> ptr;

        TreeNode::ptr               next;
        std::vector<TreeNode::ptr>  children;

        TreeNode(std::string filename, Token::ptr p) : next(NULL),filename(filename),content(p) 
        { this->children.push_back(TreeNode::ptr(NULL)); }

        void print_tree(unsigned int level=0)
        {
          if (this->content == 0) return;
          std::cout << "token type (" << this->content->TokenName[this->content->type] << ")  \t  :     "
                    << std::string(level, '\t') << this->content->value << std::endl;

          switch (this->content->type) {
            case Token::IfThen:
              {
                int i=0;
                for (auto it=this->children.begin(); it!=this->children.end(); ++it,++i) {
                  std::cout << "IF BRANCH " << i << ":" << std::endl;
                  if (*it) (*it)->print_tree(level+1);
                }
                std::cout << "END IF:" << std::endl;
              }
              break;
            default:
              break;
          }
          if (this->next) this->next->print_tree(level);
        }

        void traverse_tree(std::map<std::string,int>& defined_symbols,
                           std::set<std::string>& dependencies)
        {
          if (this->content == 0) return;

          switch (this->content->type) {
            case Token::Define:
              defined_symbols[this->content->value] = 1;
              break;
            case Token::Undefine:
              {
                auto it = defined_symbols.find(this->content->value);
                if (it != defined_symbols.end())  defined_symbols.erase(it);
              }
              break;
            case Token::Usedep:
            case Token::Incldep:
              dependencies.insert(this->content->value);
              break;
            case Token::IfThen:
              {
                // first, find out which branch to evaluate:
                int ibranch=-1, i=0;
                for (auto it=this->children.begin(); it!=this->children.end(); ++it,++i) {
                  int eval_branch = 0;
                  if (i==0)  
                    eval_branch = this->content->condition->evaluate(defined_symbols);
                  else if (!(*it)->content->condition)
                    eval_branch = 1;
                  else
                    eval_branch = (*it)->content->condition->evaluate(defined_symbols);
                  if (eval_branch) { ibranch=i; 
                     break; }
                }
                if ((ibranch >= 0) && this->children[ibranch]) 
                  this->children[ibranch]->traverse_tree(defined_symbols, dependencies);
              }
              break;
            default:
              break;
          }
          if (this->next) this->next->traverse_tree(defined_symbols, dependencies);
        }

        std::string    filename;  //< file which contains this token
        Token::ptr     content;   //< syntax token (If/Else/Use/Module/...)
    };


    FileTokenizer(const std::string& input)
      : Tokenizer(input), ilevel(0) { 
        
#line 1070 "usedep_parser.cpp"
	{
	cs = TokenParser_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 408 "usedep_parser.rl"
      }


    Token::ptr
    build_tree(std::string filename, std::map<std::string, TreeNode::ptr>& code_part, 
               TreeNode::ptr* node)
    {
       TreeNode::ptr* root = NULL;
       Token::ptr p;
       do {
         p = this->next();
         if (p->type != Token::End) {
           switch (p->type) {
             case Token::IfThen:
               {
                 // if condition: recursion for sub-branches
                 node->reset(new TreeNode(filename,p));
                 auto p1 = this->build_tree(filename,code_part, &(*node)->children[0]);
                 while ((p1->type == Token::Else) || (p1->type == Token::ElseIf))
                   {
                     (*node)->children.push_back(TreeNode::ptr(new TreeNode(filename,p1)));
                     p1 = this->build_tree(filename,code_part, &(*node)->children.back()->next);
                   }
                 node = &(*node)->next;
               }
               break;
             case Token::EndIf:
             case Token::Else:
             case Token::ElseIf:
               return p;
             case Token::Module:
               {
                 root = node; // store tree-node of the whole file
                 code_part[p->value] = TreeNode::ptr(new TreeNode(filename,p));
                 node = &code_part[p->value];
                 TreeNode::ptr src_node(code_part[filename]);
                 while (src_node) {  node->reset(new TreeNode(filename,src_node->content)); 
                                     node = &(*node)->next; src_node = src_node->next; }
               }
               break;
             case Token::EndModule:
               node = root;   // return to tree-node of the whole file
               break;
             default:
               node->reset(new TreeNode(filename,p));
               node = &(*node)->next;
           }
         }
       } while (p->type != Token::End);
    
       return p;
    }
    
private:
    int ilevel;                 //< state of ifdef nesting

    // store a token's contents (trimming whitespace):
    void capture_token(Token::ptr& token, const Token::TokenType t) 
    { transform(cbuffer.begin(),cbuffer.end(),cbuffer.begin(),tolower);
      token->value = cbuffer;  token->type = t;  }

    Token::ptr next() {
        Token::ptr token(new Token());
        token->type = Token::None;
    
        do {
            if (cs >= TokenParser_first_final)
                token->type = Token::End;
    
            
#line 1149 "usedep_parser.cpp"
	{
	int _klen;
	unsigned int _trans;
	short _widec;
	const char *_acts;
	unsigned int _nacts;
	const short *_keys;

	if ( p == pe )
		goto _test_eof;
	if ( cs == 0 )
		goto _out;
_resume:
	_acts = _TokenParser_actions + _TokenParser_from_state_actions[cs];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 ) {
		switch ( *_acts++ ) {
	case 5:
#line 1 "NONE"
	{ts = p;}
	break;
#line 1171 "usedep_parser.cpp"
		}
	}

	_widec = (*p);
	_klen = _TokenParser_cond_lengths[cs];
	_keys = _TokenParser_cond_keys + (_TokenParser_cond_offsets[cs]*2);
	if ( _klen > 0 ) {
		const short *_lower = _keys;
		const short *_mid;
		const short *_upper = _keys + (_klen<<1) - 2;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + (((_upper-_lower) >> 1) & ~1);
			if ( _widec < _mid[0] )
				_upper = _mid - 2;
			else if ( _widec > _mid[1] )
				_lower = _mid + 2;
			else {
				switch ( _TokenParser_cond_spaces[_TokenParser_cond_offsets[cs] + ((_mid - _keys)>>1)] ) {
	case 0: {
		_widec = (short)(128 + ((*p) - -128));
		if ( 
#line 252 "usedep_parser.rl"
 ilevel > 0  ) _widec += 256;
		break;
	}
				}
				break;
			}
		}
	}

	_keys = _TokenParser_trans_keys + _TokenParser_key_offsets[cs];
	_trans = _TokenParser_index_offsets[cs];

	_klen = _TokenParser_single_lengths[cs];
	if ( _klen > 0 ) {
		const short *_lower = _keys;
		const short *_mid;
		const short *_upper = _keys + _klen - 1;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + ((_upper-_lower) >> 1);
			if ( _widec < *_mid )
				_upper = _mid - 1;
			else if ( _widec > *_mid )
				_lower = _mid + 1;
			else {
				_trans += (unsigned int)(_mid - _keys);
				goto _match;
			}
		}
		_keys += _klen;
		_trans += _klen;
	}

	_klen = _TokenParser_range_lengths[cs];
	if ( _klen > 0 ) {
		const short *_lower = _keys;
		const short *_mid;
		const short *_upper = _keys + (_klen<<1) - 2;
		while (1) {
			if ( _upper < _lower )
				break;

			_mid = _lower + (((_upper-_lower) >> 1) & ~1);
			if ( _widec < _mid[0] )
				_upper = _mid - 2;
			else if ( _widec > _mid[1] )
				_lower = _mid + 2;
			else {
				_trans += (unsigned int)((_mid - _keys)>>1);
				goto _match;
			}
		}
		_trans += _klen;
	}

_match:
	_trans = _TokenParser_indicies[_trans];
_eof_trans:
	cs = _TokenParser_trans_targs[_trans];

	if ( _TokenParser_trans_actions[_trans] == 0 )
		goto _again;

	_acts = _TokenParser_actions + _TokenParser_trans_actions[_trans];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 )
	{
		switch ( *_acts++ )
		{
	case 0:
#line 240 "usedep_parser.rl"
	{ token->condition.reset(new ExpressionTokenizer(cbuffer)); }
	break;
	case 1:
#line 241 "usedep_parser.rl"
	{ cbuffer = "";  }
	break;
	case 2:
#line 242 "usedep_parser.rl"
	{ if ((*p)!='\n') cbuffer += (*p); }
	break;
	case 6:
#line 1 "NONE"
	{te = p+1;}
	break;
	case 7:
#line 278 "usedep_parser.rl"
	{act = 1;}
	break;
	case 8:
#line 279 "usedep_parser.rl"
	{act = 2;}
	break;
	case 9:
#line 280 "usedep_parser.rl"
	{act = 3;}
	break;
	case 10:
#line 282 "usedep_parser.rl"
	{act = 4;}
	break;
	case 11:
#line 286 "usedep_parser.rl"
	{act = 7;}
	break;
	case 12:
#line 294 "usedep_parser.rl"
	{act = 13;}
	break;
	case 13:
#line 283 "usedep_parser.rl"
	{te = p+1;{ if (ilevel == 0)  throw(ParserException("No preceding #if statement."));
                     capture_token(token, Token::Else);             {p++; goto _out; } }}
	break;
	case 14:
#line 285 "usedep_parser.rl"
	{te = p+1;{ ilevel--; capture_token(token, Token::EndIf);  {p++; goto _out; } }}
	break;
	case 15:
#line 286 "usedep_parser.rl"
	{te = p+1;{ capture_token(token, Token::Usedep);           {p++; goto _out; } }}
	break;
	case 16:
#line 287 "usedep_parser.rl"
	{te = p+1;{ capture_token(token, Token::Incldep);          {p++; goto _out; } }}
	break;
	case 17:
#line 290 "usedep_parser.rl"
	{te = p+1;{ capture_token(token, Token::Module);           {p++; goto _out; } }}
	break;
	case 18:
#line 291 "usedep_parser.rl"
	{te = p+1;{ capture_token(token, Token::EndModule);        {p++; goto _out; } }}
	break;
	case 19:
#line 294 "usedep_parser.rl"
	{te = p+1;{ capture_token(token, Token::Code); token->value=std::string(ts,te-ts-1);{p++; goto _out; } }}
	break;
	case 20:
#line 278 "usedep_parser.rl"
	{te = p;p--;}
	break;
	case 21:
#line 279 "usedep_parser.rl"
	{te = p;p--;{ capture_token(token, Token::IfThen); ilevel++; {p++; goto _out; } }}
	break;
	case 22:
#line 280 "usedep_parser.rl"
	{te = p;p--;{ token->condition->negate(); cbuffer = "!( " + cbuffer + " )";
                     capture_token(token, Token::IfThen);ilevel++;  {p++; goto _out; } }}
	break;
	case 23:
#line 282 "usedep_parser.rl"
	{te = p;p--;{ capture_token(token, Token::ElseIf);           {p++; goto _out; } }}
	break;
	case 24:
#line 286 "usedep_parser.rl"
	{te = p;p--;{ capture_token(token, Token::Usedep);           {p++; goto _out; } }}
	break;
	case 25:
#line 288 "usedep_parser.rl"
	{te = p;p--;{ capture_token(token, Token::Define);           {p++; goto _out; } }}
	break;
	case 26:
#line 289 "usedep_parser.rl"
	{te = p;p--;{ capture_token(token, Token::Undefine);         {p++; goto _out; } }}
	break;
	case 27:
#line 294 "usedep_parser.rl"
	{te = p;p--;{ capture_token(token, Token::Code); token->value=std::string(ts,te-ts-1);{p++; goto _out; } }}
	break;
	case 28:
#line 295 "usedep_parser.rl"
	{te = p;p--;}
	break;
	case 29:
#line 278 "usedep_parser.rl"
	{{p = ((te))-1;}}
	break;
	case 30:
#line 288 "usedep_parser.rl"
	{{p = ((te))-1;}{ capture_token(token, Token::Define);           {p++; goto _out; } }}
	break;
	case 31:
#line 1 "NONE"
	{	switch( act ) {
	case 0:
	{{cs = 0;goto _again;}}
	break;
	case 2:
	{{p = ((te))-1;} capture_token(token, Token::IfThen); ilevel++; {p++; goto _out; } }
	break;
	case 3:
	{{p = ((te))-1;} token->condition->negate(); cbuffer = "!( " + cbuffer + " )";
                     capture_token(token, Token::IfThen);ilevel++;  {p++; goto _out; } }
	break;
	case 4:
	{{p = ((te))-1;} capture_token(token, Token::ElseIf);           {p++; goto _out; } }
	break;
	case 7:
	{{p = ((te))-1;} capture_token(token, Token::Usedep);           {p++; goto _out; } }
	break;
	case 13:
	{{p = ((te))-1;} capture_token(token, Token::Code); token->value=std::string(ts,te-ts-1);{p++; goto _out; } }
	break;
	default:
	{{p = ((te))-1;}}
	break;
	}
	}
	break;
#line 1410 "usedep_parser.cpp"
		}
	}

_again:
	_acts = _TokenParser_actions + _TokenParser_to_state_actions[cs];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 ) {
		switch ( *_acts++ ) {
	case 3:
#line 1 "NONE"
	{ts = 0;}
	break;
	case 4:
#line 1 "NONE"
	{act = 0;}
	break;
#line 1427 "usedep_parser.cpp"
		}
	}

	if ( cs == 0 )
		goto _out;
	if ( ++p != pe )
		goto _resume;
	_test_eof: {}
	if ( p == eof )
	{
	if ( _TokenParser_eof_trans[cs] > 0 ) {
		_trans = _TokenParser_eof_trans[cs] - 1;
		goto _eof_trans;
	}
	}

	_out: {}
	}

#line 478 "usedep_parser.rl"
            
            if (cs == TokenParser_error) 
              throw(ParserException("Tokenizer parsing failed."));

        } while (token->type == Token::None);
    
        // consistency checks
        if ((token->type == Token::End) && (this->ilevel != 0))
          throw(ParserException("Unbalanced ifdef-end construct."));
        return token;
    }
};


// --- Auxiliary routine: Check if file exists
inline bool fileExists(const std::string& filename) {
  struct stat st_buf;
  int status = stat(filename.c_str(), &st_buf);
  return (status == 0) && (S_ISREG (st_buf.st_mode)); 
}
// --- Auxiliary routine: Check if string ends with a substring
inline bool ends_with(std::string const& value, std::string const& ending)
{
  if (ending.size() > value.size()) return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}
// --- Auxiliary routine: strip file name from path
inline std::string wout_path(std::string const& filename)
{
  auto pos = filename.find_last_of("/");
  if (pos == std::string::npos)  
    return filename;
  else 
    return filename.substr(pos+1,filename.length());
}
// --- Auxiliary routine: transform file extension to ".o"
inline std::string obj_name(std::string const& filename, std::string const& objprefix)
{
  auto pos = filename.find_last_of(".");
  if (pos == std::string::npos)  
    return objprefix+filename;
  else 
    return objprefix+filename.substr(0,pos)+".o";
}


// --------------------------------------------------------------------------------
// Main routine
// --------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  if (argc<4) {
    std::cout << "Usage: " << argv[0] 
              << "<root> <defined_symbols> <filenames> [-v] [-m] [-objprefix=xxx]" << std::endl; 
     return -1;
  }
  std::string search_root(argv[1]), defined_symbols_str(argv[2]), filelist_str(argv[3]);
  // optional command-line arguments:
  std::map<std::string, std::string> optional_args;
  for (int i=4; i<argc; ++i) {
    std::string arg = std::string(argv[i]);
    if (arg[0] == '-')  arg.erase(0,1);
    std::string value = "";
    auto pos = arg.find("=");
    if (pos != std::string::npos) { value=arg.substr(pos+1,arg.length());arg=arg.substr(0,pos); }
    optional_args[arg] = value;
  }

  bool verbose           = (optional_args.find("v") != optional_args.end());
  bool print_convex_hull = (optional_args.find("m") == optional_args.end());
  std::string objprefix = "";
  if (optional_args.find("objprefix") != optional_args.end()) objprefix = optional_args["objprefix"];

  // parse lists of defined symbols and input files from command line
  typedef std::istream_iterator<std::string> str_iterator;
  std::istringstream           defined_symbols_sstr(defined_symbols_str);
  std::vector<std::string>     tokens{str_iterator{defined_symbols_sstr}, str_iterator{}};
  std::istringstream           filelist_sstr(filelist_str);
  std::vector<std::string>     files{str_iterator{filelist_sstr}, str_iterator{}};

  std::map< std::string, int > defined_symbols;
  for (auto it=tokens.begin(); it!=tokens.end(); ++it) defined_symbols[*it] = 1;

  std::map< std::string, FileTokenizer::TreeNode::ptr > code_part;

  // process all files; build a  map of code parts:  MODULE NAME -> TOKEN TREE
  for (auto file_it=files.begin(); file_it!=files.end(); ++file_it) 
  {
    if (!fileExists(*file_it)) {
      std::cout << "Error! File \"" << *file_it << "\" does not exist!" << std::endl;  return -1;
    }
    // read the whole file into a string
    std::stringstream buffer;
    buffer << std::ifstream(*file_it).rdbuf();

    FileTokenizer tokenizer(buffer.str());
    try {
      if (verbose)  std::cout << "processing \"" << *file_it << "\"" << std::endl;
      code_part[*file_it] = FileTokenizer::TreeNode::ptr(new FileTokenizer::TreeNode(*file_it,NULL));
      tokenizer.build_tree(*file_it, code_part, &code_part[*file_it]); 
    } 
    catch (std::exception& e) {
      // provide detailed context for the error
      std::cout << std::endl << "Error (" << *file_it << "): " << e.what() << std::endl;
      if (tokenizer.ts != NULL) {
        const char *line_start = tokenizer.ts, *line_end = tokenizer.ts;
        while ((line_start > tokenizer.buffer.c_str()) && (*line_start != '\n'))  --line_start;
        while ((line_end   < tokenizer.pe)             && (*line_end   != '\n'))  ++line_end;
        if (line_end > line_start)
          std::cout << "Error somewhere around " << std::endl
                    << "\"" << std::string(++line_start,line_end) << " [...]\"" << std::endl;
      }
    }
  }

  // traverse token tree and find dependencies of "search_root" file
  std::set<std::string> visited;
  std::vector<std::string> list; list.push_back(search_root);
  std::map< std::string,std::set<std::string> > make_dependencies;
  std::map< std::string, std::string>          include_file;
  while (!list.empty()) {
    auto item = list.back();  list.pop_back();
    auto cp = code_part.find(item);  auto cp0=cp;
    if (cp == code_part.end()) { // alternatively prepend symbol by base path
      for (auto it2=files.begin(); it2!=files.end(); ++it2)
        if (ends_with(*it2, "/"+item)) { 
          include_file[item] = *it2;
          item=*it2; visited.insert(*it2); cp=code_part.find(item); break;
        }
    }
    if (cp != code_part.end()) {
      std::set<std::string> dependencies;
      cp->second->traverse_tree(defined_symbols, dependencies);
      // first, construct a list of all dependencies ("convex hull"):
      for (auto it2=dependencies.begin(); it2!=dependencies.end(); ++it2)
        if (visited.find(*it2) == visited.end()) 
          { list.push_back(*it2); visited.insert(*it2); }

      // besides, for each module (not "include") build a list of direct dependencies:
      if (cp0 != code_part.end()) {
        const std::string& fname = cp->second->filename;
        auto make_dep = make_dependencies.find(fname);
        if (make_dep == make_dependencies.end())  
          { make_dependencies[fname] = std::set<std::string>(); make_dep = make_dependencies.find(fname); }
        make_dep->second.insert(dependencies.begin(), dependencies.end());
      }
    }
    else if (verbose)  
      std::cout << "USE/include outside of scope: \"" << item << "\"" << std::endl;
  }

  if (print_convex_hull) 
    {
      // output variant I: print out files of union ("convex hull") 
      //                   of all dependencies:
      if (verbose)  std::cout << std::endl << "RESULT: \"" << search_root << "\" depends on the following files" 
                              << std::endl << "        if \"" << defined_symbols_str << "\" is set:" << std::endl;
      std::set< std::string > all_deps;
      for (auto it=visited.begin(); it!=visited.end(); ++it)
        if (code_part.find(*it) != code_part.end())  all_deps.insert(code_part[*it]->filename);
      for (auto it=all_deps.begin(); it!=all_deps.end(); ++it)
        std::cout << *it << std::endl;
    }
  else
    {
      // output variant II: print out immediate dependencies for each
      //                    element of convex hull
      for (auto it=make_dependencies.begin(); it!=make_dependencies.end(); ++it) {
        if (code_part.find(it->first) != code_part.end()) {

          std::set<std::string> local_list;
          for (auto it2=it->second.begin(); it2!=it->second.end(); ++it2) 
            {
              std::string name = (*it2);
              auto cp_name = code_part.find(name);
              if (cp_name != code_part.end()) 
                // direct dependencies which are no includes are
                // resolved into object names:
                local_list.insert(obj_name(wout_path(cp_name->second->filename), objprefix));
              else
                // while includes are prepended by the full path name:
                {
                  auto incl_name = include_file.find(name);
                  if (incl_name != include_file.end())  local_list.insert(incl_name->second);
                }
            }

          // print out the direct dependencies:
          std::cout << obj_name(wout_path(it->first),objprefix) << " : " << it->first;
          for (auto it2=local_list.begin(); it2!=local_list.end(); ++it2)
            std::cout << " \\" << std::endl << *it2;
          std::cout << std::endl << std::endl;
        }
      }
    }

  return 0;
}
