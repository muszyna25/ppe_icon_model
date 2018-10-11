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

%%{
  machine ExpressionParser;

    action clear   { cbuffer = "";  }
    action append  { if (fc!='\n') cbuffer += fc; }

    # store operator 
    action store_op { 
      /*  While the stack is not empty and an operator is at the top
          and the operator at the top is higher priority than the item
          then pop the operator on the top of the stack, add the
          popped operator to the queue. */
      int pr = priority(cbuffer);
      while ((!ostack.empty()) && (ostack.back().type==OPERATOR) &&
	     (pr < priority(ostack.back().val)))
        { queue.push_back(ostack.back()); ostack.pop_back(); }
      ostack.push_back( t_item(OPERATOR, cbuffer) );
     }
    # store parentheses
    action store_paren {
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
    }
  
  parenthesis    = ('('|')') >clear @append ;
  value          = (digit | '-' digit) >clear $append  digit*         @append ;
  varname        = (alpha | [_])       >clear @append  (alnum | [_])* @append ;
  defined_token  = /DEFINED/i space* ;
  operator_token = ( '&' | '|' | '>' | '<' | '!' )+ >clear @append ;
  
  main := |*
    parenthesis                                 => store_paren;
    value                                       => { queue.push_back( t_item(VALUE,    cbuffer) ); };
    varname                                     => { queue.push_back( t_item(VARIABLE, cbuffer) ); };
    defined_token '(' space* varname space* ')' => { queue.push_back( t_item(VARIABLE, cbuffer) ); };
    defined_token space* varname space*         => { queue.push_back( t_item(VARIABLE, cbuffer) ); };
    operator_token                              => store_op;
    space+;
  *| ;

}%%

%% write data nofinal;


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
    %%write init;    
    %%write exec;

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

%%{
  machine TokenParser;

  action store_condition   { token->condition.reset(new ExpressionTokenizer(cbuffer)); }
  action clear             { cbuffer = "";  }
  action append            { if (fc!='\n') cbuffer += fc; }

  main := |*
    # Preprocessor ifdef pattern
    # nested ifdef's are handled via semantic conditions:
    hash        = '#' ' '* ;
    condition   = (alnum | [()&\-| _><!])+ >clear @append %store_condition ;
    open_ifdef  = (hash 'ifdef' | '#if' space) space* condition ;
    open_ifndef = hash 'ifndef' space* condition ;
    open_elif   = hash 'elif'   space* condition ;
    close_ifdef = ( hash 'endif' when { ilevel > 0 } ) >clear ;
    else_ifdef  = hash 'else' >clear ;

    # MODULE/END MODULE name pattern
    valid_name  = (alpha|[_]) >clear @append (alnum | [_])* @append ;
    module      = space* /MODULE/i space+ valid_name [ ]* '\n' ;
    endmodule   = space* /END/i [ \t]* /MODULE/i space+ valid_name [ ]* '\n' ;

    # USE dependency pattern
    comment     = '!' [^\n]* ; 
    use_spec    = (',' space* [^\n]*) ;
    usedep      = space* /USE/i space+ valid_name space* (comment|use_spec)? '\n' ;

    # INCLUDE dependency pattern
    filename    = alpha >clear @append (alnum | [._])* @append ;
    includedep  = hash 'include' space* ['"] filename ['"] ;

    # DEFINE/UNDEFINE pattern
    define       = hash 'define' space+ valid_name ;
    undefine     = hash 'undef'  space+ valid_name ;
    cmplx_define = hash 'define' space+ valid_name'(' alnum+ ')' space+ hash  ;

    # Code pattern (the rest)
    string_literal = (('"' [^"]* '"') |  ("'" [^']* "'")) ;
    code           = (([^\n#!'"]) | string_literal )* >clear comment?  '\n' ;

    space+;
    open_ifdef  => { capture_token(token, Token::IfThen); ilevel++; fbreak; };
    open_ifndef => { token->condition->negate(); cbuffer = "!( " + cbuffer + " )";
                     capture_token(token, Token::IfThen);ilevel++;  fbreak; };
    open_elif   => { capture_token(token, Token::ElseIf);           fbreak; };
    else_ifdef  => { if (ilevel == 0)  throw(ParserException("No preceding #if statement."));
                     capture_token(token, Token::Else);             fbreak; };
    close_ifdef => { ilevel--; capture_token(token, Token::EndIf);  fbreak; };
    usedep      => { capture_token(token, Token::Usedep);           fbreak; };
    includedep  => { capture_token(token, Token::Incldep);          fbreak; };
    define      => { capture_token(token, Token::Define);           fbreak; };
    undefine    => { capture_token(token, Token::Undefine);         fbreak; };
    module      => { capture_token(token, Token::Module);           fbreak; };
    endmodule   => { capture_token(token, Token::EndModule);        fbreak; }; 

    # note: ordering the code token ("the rest") last is important:
    code        => { capture_token(token, Token::Code); token->value=std::string(ts,te-ts-1);fbreak; };
    cmplx_define; # this captures the complex define "#define STRINGIFY(a) #a"
  *|;
}%%

%% write data;


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
        %%write init;
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
    
            %%write exec;
            
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
