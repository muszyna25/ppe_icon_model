
#line 1 "xml_scan.rl"
/* -------------------------------------------------------------------
 * Generic light-weight XML parser with some access functions for grid
 * description; uses Ragel finite state machine compiler.
 *
 * 2013-04-05 : F. Prill, DWD
 * -------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>


/*
 * Single attribute of an XML tag.
 */
struct t_attr {
  char              *key, *value;  /* key/value pair of attribute */
  struct t_attr     *next;         /* next attribute */
};

/*
 * Single XML tag.
 */
struct t_element {
  char              *name, *body;  /* element name, body text */
  int               level;         /* level in XML tree */
  struct t_attr     *attr;         /* first attribute */
  struct t_element  *child, *next; /* first child element; sibling element */
};

struct t_element *xmltree = NULL; /* head of XML tree */


/* --------------------------------------------------------------
 * Simple XML parser
 * -------------------------------------------------------------- */


#line 111 "xml_scan.rl"



#line 48 "xml_scan.c"
static const char _xmlparse_actions[] = {
	0, 1, 0, 1, 1, 1, 2, 1, 
	3, 1, 4, 1, 5, 1, 6, 1, 
	7, 2, 1, 4, 2, 2, 4
};

static const char _xmlparse_key_offsets[] = {
	0, 0, 4, 9, 14, 20, 26, 30, 
	35, 36, 41, 45, 49, 50, 54, 55, 
	60, 65, 71, 77, 81, 86, 87, 92, 
	96, 100, 101, 105, 106, 107, 108, 111
};

static const char _xmlparse_trans_keys[] = {
	32, 60, 9, 13, 32, 47, 62, 9, 
	13, 32, 47, 62, 9, 13, 32, 47, 
	61, 62, 9, 13, 32, 47, 61, 62, 
	9, 13, 32, 61, 9, 13, 32, 34, 
	39, 9, 13, 34, 32, 47, 62, 9, 
	13, 32, 62, 9, 13, 32, 62, 9, 
	13, 39, 32, 60, 9, 13, 60, 32, 
	47, 62, 9, 13, 32, 47, 62, 9, 
	13, 32, 47, 61, 62, 9, 13, 32, 
	47, 61, 62, 9, 13, 32, 61, 9, 
	13, 32, 34, 39, 9, 13, 34, 32, 
	47, 62, 9, 13, 32, 62, 9, 13, 
	32, 62, 9, 13, 60, 32, 47, 9, 
	13, 62, 62, 39, 32, 9, 13, 0
};

static const char _xmlparse_single_lengths[] = {
	0, 2, 3, 3, 4, 4, 2, 3, 
	1, 3, 2, 2, 1, 2, 1, 3, 
	3, 4, 4, 2, 3, 1, 3, 2, 
	2, 1, 2, 1, 1, 1, 1, 0
};

static const char _xmlparse_range_lengths[] = {
	0, 1, 1, 1, 1, 1, 1, 1, 
	0, 1, 1, 1, 0, 1, 0, 1, 
	1, 1, 1, 1, 1, 0, 1, 1, 
	1, 0, 1, 0, 0, 0, 1, 0
};

static const unsigned char _xmlparse_index_offsets[] = {
	0, 0, 4, 9, 14, 20, 26, 30, 
	35, 37, 42, 46, 50, 52, 56, 58, 
	63, 68, 74, 80, 84, 89, 91, 96, 
	100, 104, 106, 110, 112, 114, 116, 119
};

static const char _xmlparse_indicies[] = {
	0, 2, 0, 1, 2, 1, 1, 2, 
	3, 4, 5, 6, 4, 3, 8, 9, 
	1, 10, 8, 7, 11, 1, 12, 1, 
	11, 7, 13, 14, 13, 1, 14, 15, 
	16, 14, 1, 18, 17, 8, 9, 10, 
	8, 1, 19, 20, 19, 1, 21, 10, 
	21, 1, 18, 22, 24, 25, 24, 23, 
	26, 23, 25, 28, 1, 25, 27, 29, 
	30, 31, 29, 27, 33, 34, 1, 35, 
	33, 32, 36, 1, 37, 1, 36, 32, 
	38, 39, 38, 1, 39, 40, 41, 39, 
	1, 43, 42, 33, 34, 35, 33, 1, 
	44, 45, 44, 1, 46, 35, 46, 1, 
	47, 1, 47, 28, 47, 1, 1, 48, 
	49, 48, 43, 50, 51, 51, 1, 1, 
	0
};

static const char _xmlparse_trans_targs[] = {
	1, 0, 2, 3, 4, 10, 30, 5, 
	4, 10, 30, 6, 7, 6, 7, 8, 
	12, 8, 9, 11, 30, 11, 12, 14, 
	13, 15, 15, 16, 27, 17, 23, 25, 
	18, 17, 23, 25, 19, 20, 19, 20, 
	21, 29, 21, 22, 24, 25, 24, 26, 
	28, 31, 29, 30
};

static const char _xmlparse_trans_actions[] = {
	0, 0, 0, 1, 3, 3, 17, 1, 
	0, 0, 9, 11, 11, 0, 0, 0, 
	0, 1, 13, 5, 20, 0, 1, 1, 
	0, 0, 15, 1, 0, 3, 3, 17, 
	1, 0, 0, 9, 11, 11, 0, 0, 
	0, 0, 1, 13, 5, 20, 0, 0, 
	0, 7, 1, 0
};

static const int xmlparse_start = 1;
static const int xmlparse_error = 0;

static const int xmlparse_en_elementBody = 13;
static const int xmlparse_en_main = 1;


#line 114 "xml_scan.rl"

#define BUFSIZE 128

/*
 * Opens an XML file and parses its contents into a tree-like data
 * structure.
 */
void parseXML(const char* filename)
{
  static char buf[BUFSIZE];
  static char buf2[BUFSIZE];
  int   cs, have  = 0;
  char *ts, *te   = 0;
  int   done      = 0;
  FILE *file      = fopen(filename , "r");
  int   elevel    = 0;
  struct t_element* tnode = NULL;

  
#line 166 "xml_scan.c"
	{
	cs = xmlparse_start;
	}

#line 133 "xml_scan.rl"

  int s = 0;
  while ( !done ) {
    char *p = buf + have, *pe, *eof = 0;
    int len, space = BUFSIZE - have;
    if ( space == 0 ) {
      fprintf(stderr, "OUT OF BUFFER SPACE\n" );
      exit(1);
    }
    len = fread( p, 1, space, file );
    pe = p + len;
    if ( len < space ) {
      eof  = pe;
      done = 1;
    }
    
    
#line 189 "xml_scan.c"
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
	_keys = _xmlparse_trans_keys + _xmlparse_key_offsets[cs];
	_trans = _xmlparse_index_offsets[cs];

	_klen = _xmlparse_single_lengths[cs];
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

	_klen = _xmlparse_range_lengths[cs];
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
	_trans = _xmlparse_indicies[_trans];
	cs = _xmlparse_trans_targs[_trans];

	if ( _xmlparse_trans_actions[_trans] == 0 )
		goto _again;

	_acts = _xmlparse_actions + _xmlparse_trans_actions[_trans];
	_nacts = (unsigned int) *_acts++;
	while ( _nacts-- > 0 )
	{
		switch ( *_acts++ )
		{
	case 0:
#line 43 "xml_scan.rl"
	{
	  buf2[s] = *p; s++;
	}
	break;
	case 1:
#line 46 "xml_scan.rl"
	{
	  elevel++; 
	  if (tnode == NULL) {
	    xmltree = malloc(sizeof(struct t_element));
	    tnode   = xmltree;
	  } 
	  else  {
	    if (elevel == tnode->level) {
	      tnode->next = malloc(sizeof(struct t_element));
	      tnode = tnode->next;
	    }
	    else {
	      tnode->child = malloc(sizeof(struct t_element));
	      tnode = tnode->child;
	    }
	  }
          tnode->level = elevel;
	  tnode->child = NULL;
	  tnode->next  = NULL;
	  tnode->attr  = NULL;

	  tnode->name = malloc(sizeof(char)*(s+1));
	  strncpy(tnode->name, buf2, s);
	  s = 0;
	}
	break;
	case 2:
#line 71 "xml_scan.rl"
	{ elevel--; {cs = 13; goto _again;} }
	break;
	case 3:
#line 72 "xml_scan.rl"
	{ elevel--; {cs = 13; goto _again;} }
	break;
	case 4:
#line 73 "xml_scan.rl"
	{ {cs = 13; goto _again;} }
	break;
	case 5:
#line 74 "xml_scan.rl"
	{
	  struct t_attr *cur_attr = tnode->attr;
	  if (cur_attr == NULL) {
	    tnode->attr = malloc(sizeof(struct t_attr));
	    cur_attr = tnode->attr;
	  }
	  else {
	    while (cur_attr->next != NULL) cur_attr = cur_attr->next;
	    cur_attr->next = malloc(sizeof(struct t_attr));
	    cur_attr = cur_attr->next;
	  }
	  cur_attr->value = NULL;
	  cur_attr->next  = NULL;
	  cur_attr->key   = malloc(sizeof(char)*(s+1));
	  strncpy(cur_attr->key, buf2, s);
	  s = 0;
	}
	break;
	case 6:
#line 91 "xml_scan.rl"
	{
	  struct t_attr *cur_attr = tnode->attr;
	  while (cur_attr->next != NULL) cur_attr = cur_attr->next;
	  cur_attr->value = malloc(sizeof(char)*(s+1));
	  strncpy(cur_attr->value, buf2, s);
	  s = 0;
	}
	break;
	case 7:
#line 98 "xml_scan.rl"
	{
	  tnode->body = malloc(sizeof(char)*(s+1));
	  strncpy(tnode->body, buf2, s);
	  s = 0;
	}
	break;
#line 347 "xml_scan.c"
		}
	}

_again:
	if ( cs == 0 )
		goto _out;
	if ( ++p != pe )
		goto _resume;
	_test_eof: {}
	_out: {}
	}

#line 150 "xml_scan.rl"

    if ( cs == xmlparse_error ) {
      fprintf(stderr, "PARSE ERROR\n" );
      break;
    }
    if ( ts == 0 ) 
      have = 0;
    else {
      /* There is a prefix to preserve, shift it over.  */
      have = pe - ts;
      memmove( buf, ts, have );
      te = buf + (te-ts);
      ts = buf;
    }
  }
  fclose(file);
}


/* --------------------------------------------------------------
 * Generic routines for accessing parsed XML tree
 * -------------------------------------------------------------- */

void print_xml_tree(struct t_element *node) {
  struct t_element *tnode = NULL;
  struct t_attr    *tattr = NULL;
  
  if (node->body != NULL)   
    printf("element %d/%s (%s)\n", node->level, node->name, node->body);
  else
    printf("element %d/%s\n", node->level, node->name);

  printf("attributes:\n");
  tattr = node->attr;
  while (tattr != NULL) {
    printf("  %s: %s\n", tattr->key, tattr->value);
    tattr = tattr->next;
  }

  printf("siblings:\n");
  tnode = node->next;
  while (tnode != NULL) {
    print_xml_tree(tnode);
    tnode = tnode->next;
  }
  printf("children:\n");
  tnode = node->child;
  while (tnode != NULL) {
    print_xml_tree(tnode);
    tnode = tnode->child;
  }
}


/*
 * Find first XML tag in tree which matches a given name and
 * (optional) different criteria checked by a given comparison
 * function.
 */
struct t_element* find_first_element(struct t_element *node, 
				     const  char      *name,
				     int (*cmp_function)(struct t_element*, void*),
				     void             *context) {
  struct t_element *tnode = NULL, *result = NULL;

  if (strcmp(node->name, name) == 0) {
    if (cmp_function != NULL) {
      if (cmp_function(node, context)) return node;
    }
    else 
      return node;
  }
  /* traverse siblings */
  tnode = node->next;
  while (tnode != NULL) {
    result = find_first_element(tnode, name, cmp_function,context);
    if (result != NULL) return result;
    tnode = tnode->next;
  }
  /* traverse children */
  tnode = node->child;
  while (tnode != NULL) {
    result = find_first_element(tnode, name, cmp_function,context);
    if (result != NULL) return result;
    tnode = tnode->child;
  }
  return NULL;
}


/*
 * @return first XML tag attribute matching a given name.
 */
struct t_attr* find_first_attribute(struct t_element *node, const char* key) {
  struct t_attr *tnode = NULL;

  tnode = node->attr;
  while (tnode != NULL) {
    if (strcmp(tnode->key, key) == 0) return tnode;
    tnode = tnode->next;
  }
  return NULL;
}


/* --------------------------------------------------------------
 * Application of XML parser to <grid> :: <uri>, <description>
 * -------------------------------------------------------------- */

struct t_grid_attr {
  int   number_of_grid_used, centre, subcentre;
  char* type;
};


int match_grid(struct t_element *node, void *context) {
  struct t_grid_attr *ta = (struct t_grid_attr *) context;
  struct t_attr* attr = NULL;
  attr = find_first_attribute(node, "number_of_grid_used");
  if (atoi(attr->value) != ta->number_of_grid_used) return 0;
  attr = find_first_attribute(node, "centre");
  if (atoi(attr->value) != ta->centre) return 0;
  attr = find_first_attribute(node, "subcentre");
  if (atoi(attr->value) != ta->subcentre) return 0;
  attr = find_first_attribute(node, "type");
  if (strcmp(attr->value, ta->type) != 0) return 0;
  return 1;
}


const char* get_grid_uri(struct t_grid_attr *grid_attr) {
  /* step 1: find <grid> tag: */
  struct t_element *grid_node = find_first_element(xmltree, "grid", 
						   &match_grid, grid_attr);
  if (grid_node == NULL) return NULL;
  /* step 2: get corresponding URI tag: */
  struct t_element *uri_node = find_first_element(grid_node, "uri", NULL, NULL);
  assert(uri_node       != NULL);
  assert(uri_node->body != NULL);
  return uri_node->body;
}


const char* get_grid_description(struct t_grid_attr *grid_attr) {
  /* step 1: find <grid> tag: */
  struct t_element *grid_node = find_first_element(xmltree, "grid", 
						   &match_grid, grid_attr);
  if (grid_node == NULL) return NULL;
  /* step 2: get corresponding URI tag: */
  struct t_element *uri_node = find_first_element(grid_node, "description", NULL, NULL);
  assert(uri_node       != NULL);
  return uri_node->body;
}


/* --------------------------------------------------------------
 * Main program
 * -------------------------------------------------------------- */

int main()
{
  /* we just read an XML file and look for a given grid URI and
     description. */
  parseXML("test.xml");

  struct t_grid_attr grid_attr = {.number_of_grid_used = 137, 
				  .centre              =  78,
				  .subcentre           =   0,
				  .type                = "horizontal"};
  printf("URI: %s\n", get_grid_uri(&grid_attr));
  printf("     \"%s\"\n", get_grid_description(&grid_attr));
  return 0;
}

