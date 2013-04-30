
#line 1 "xml_scan.rl"
#ifndef XML_SCAN_H
#define XML_SCAN_H

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
int elnum = 0;

/* --------------------------------------------------------------
 * Simple XML parser
 * -------------------------------------------------------------- */


#line 123 "xml_scan.rl"



#line 51 "xml_scan.h"
static const char _xmlparse_actions[] = {
	0, 1, 0, 1, 1, 1, 2, 1, 
	3, 1, 4, 1, 5, 1, 6, 1, 
	7, 2, 0, 1, 2, 1, 2, 2, 
	1, 4, 2, 2, 4, 2, 5, 2, 
	3, 1, 2, 4
};

static const unsigned char _xmlparse_key_offsets[] = {
	0, 0, 4, 9, 15, 22, 28, 32, 
	37, 38, 44, 48, 52, 53, 59, 64, 
	70, 74, 75, 80, 86, 93, 99, 103, 
	108, 109, 115, 119, 123, 124, 128, 129, 
	130, 131, 137, 142, 148, 151
};

static const char _xmlparse_trans_keys[] = {
	32, 60, 9, 13, 32, 47, 62, 9, 
	13, 32, 47, 62, 63, 9, 13, 32, 
	47, 61, 62, 63, 9, 13, 32, 47, 
	61, 62, 9, 13, 32, 61, 9, 13, 
	32, 34, 39, 9, 13, 34, 32, 47, 
	62, 63, 9, 13, 32, 62, 9, 13, 
	32, 62, 9, 13, 39, 32, 47, 61, 
	62, 9, 13, 32, 61, 62, 9, 13, 
	32, 47, 62, 63, 9, 13, 32, 60, 
	9, 13, 60, 32, 47, 62, 9, 13, 
	32, 47, 62, 63, 9, 13, 32, 47, 
	61, 62, 63, 9, 13, 32, 47, 61, 
	62, 9, 13, 32, 61, 9, 13, 32, 
	34, 39, 9, 13, 34, 32, 47, 62, 
	63, 9, 13, 32, 62, 9, 13, 32, 
	62, 9, 13, 60, 32, 47, 9, 13, 
	62, 62, 39, 32, 47, 61, 62, 9, 
	13, 32, 61, 62, 9, 13, 32, 47, 
	62, 63, 9, 13, 32, 9, 13, 0
};

static const char _xmlparse_single_lengths[] = {
	0, 2, 3, 4, 5, 4, 2, 3, 
	1, 4, 2, 2, 1, 4, 3, 4, 
	2, 1, 3, 4, 5, 4, 2, 3, 
	1, 4, 2, 2, 1, 2, 1, 1, 
	1, 4, 3, 4, 1, 0
};

static const char _xmlparse_range_lengths[] = {
	0, 1, 1, 1, 1, 1, 1, 1, 
	0, 1, 1, 1, 0, 1, 1, 1, 
	1, 0, 1, 1, 1, 1, 1, 1, 
	0, 1, 1, 1, 0, 1, 0, 0, 
	0, 1, 1, 1, 1, 0
};

static const unsigned char _xmlparse_index_offsets[] = {
	0, 0, 4, 9, 15, 22, 28, 32, 
	37, 39, 45, 49, 53, 55, 61, 66, 
	72, 76, 78, 83, 89, 96, 102, 106, 
	111, 113, 119, 123, 127, 129, 133, 135, 
	137, 139, 145, 150, 156, 159
};

static const char _xmlparse_indicies[] = {
	0, 2, 0, 1, 2, 1, 1, 2, 
	3, 4, 5, 6, 7, 4, 3, 9, 
	10, 1, 11, 12, 9, 8, 13, 1, 
	14, 1, 13, 8, 15, 16, 15, 1, 
	16, 17, 18, 16, 1, 20, 19, 9, 
	10, 11, 10, 9, 1, 21, 22, 21, 
	1, 23, 11, 23, 1, 20, 24, 25, 
	1, 14, 22, 25, 8, 26, 16, 11, 
	26, 1, 27, 5, 28, 7, 27, 3, 
	30, 31, 30, 29, 32, 29, 31, 34, 
	1, 31, 33, 35, 36, 37, 38, 35, 
	33, 40, 41, 1, 42, 43, 40, 39, 
	44, 1, 45, 1, 44, 39, 46, 47, 
	46, 1, 47, 48, 49, 47, 1, 51, 
	50, 40, 41, 42, 41, 40, 1, 52, 
	53, 52, 1, 54, 42, 54, 1, 55, 
	1, 55, 34, 55, 1, 1, 56, 57, 
	56, 51, 58, 59, 1, 45, 53, 59, 
	39, 60, 47, 42, 60, 1, 61, 36, 
	62, 38, 61, 33, 63, 63, 1, 1, 
	0
};

static const char _xmlparse_trans_targs[] = {
	1, 0, 2, 3, 4, 10, 36, 15, 
	5, 4, 10, 36, 13, 6, 7, 6, 
	7, 8, 12, 8, 9, 11, 36, 11, 
	12, 14, 14, 4, 36, 17, 16, 18, 
	18, 19, 30, 20, 26, 28, 35, 21, 
	20, 26, 28, 33, 22, 23, 22, 23, 
	24, 32, 24, 25, 27, 28, 27, 29, 
	31, 37, 32, 34, 34, 20, 28, 36
};

static const char _xmlparse_trans_actions[] = {
	0, 0, 0, 1, 3, 3, 23, 17, 
	1, 0, 0, 9, 1, 11, 11, 0, 
	0, 0, 0, 1, 13, 5, 26, 0, 
	1, 29, 0, 20, 32, 1, 0, 0, 
	15, 1, 0, 3, 3, 23, 17, 1, 
	0, 0, 9, 1, 11, 11, 0, 0, 
	0, 0, 1, 13, 5, 26, 0, 0, 
	0, 7, 1, 29, 0, 20, 32, 0
};

static const int xmlparse_start = 1;
static const int xmlparse_error = 0;

static const int xmlparse_en_elementBody = 16;
static const int xmlparse_en_main = 1;


#line 126 "xml_scan.rl"

#define BUFSIZE       128
#define MAX_XMLLEVELS 128

/*
 * Opens an XML file and parses its contents into a tree-like data
 * structure.
 */
void parseXML(const char* filename)
{
  static char buf[BUFSIZE];
  static char buf2[BUFSIZE];
  int   cs, have  = 0;
  char *ts        = 0;
  char *te        = 0;
  int   done      = 0;
  FILE *file      = fopen(filename , "r");
  int   elevel    = 0;
  struct t_element* tnode[MAX_XMLLEVELS];
  int i;

  for (i=0; i< MAX_XMLLEVELS; i++) tnode[i] = NULL;

  
#line 192 "xml_scan.h"
	{
	cs = xmlparse_start;
	}

#line 150 "xml_scan.rl"

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
    if ((len < space) || (len == 0)) {
      eof  = pe;
      done = 1;
    }
    
    
#line 215 "xml_scan.h"
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
#line 46 "xml_scan.rl"
	{
	  buf2[s] = *p; s++;
	}
	break;
	case 1:
#line 49 "xml_scan.rl"
	{
	  elevel++; 
	  struct t_element* new_node;
	  if (tnode[elevel] == NULL) {
	    if (elevel == 1) 
	      {
		xmltree       = calloc(1, sizeof(struct t_element));
		tnode[elevel] = xmltree;
		new_node = tnode[elevel];
	      } 
	    else
	      {
		tnode[elevel-1]->child = calloc(1, sizeof(struct t_element));
		tnode[elevel] = tnode[elevel-1]->child;
		new_node = tnode[elevel];
		elnum++;
	      }
	    elnum++;
	  } 
	  else  {
	    tnode[elevel]->next = calloc(1, sizeof(struct t_element));
	    tnode[elevel] = tnode[elevel]->next;
	    new_node = tnode[elevel];
	    elnum++;
	  }
          new_node->level = elevel;
	  new_node->child = NULL;
	  new_node->next  = NULL;
	  new_node->attr  = NULL;

	  new_node->name = calloc((s+1), sizeof(char));
	  strncpy(new_node->name, buf2, s);
	  s = 0;
	}
	break;
	case 2:
#line 83 "xml_scan.rl"
	{ tnode[elevel+1] = NULL; elevel--; s = 0; {cs = 16; goto _again;} }
	break;
	case 3:
#line 84 "xml_scan.rl"
	{ tnode[elevel+1] = NULL; elevel--; s = 0; {cs = 16; goto _again;} }
	break;
	case 4:
#line 85 "xml_scan.rl"
	{ {cs = 16; goto _again;} }
	break;
	case 5:
#line 86 "xml_scan.rl"
	{
	  struct t_attr *cur_attr = tnode[elevel]->attr;
	  if (cur_attr == NULL) {
	    tnode[elevel]->attr = malloc(sizeof(struct t_attr));
	    cur_attr = tnode[elevel]->attr;
	  }
	  else {
	    while (cur_attr->next != NULL) cur_attr = cur_attr->next;
	    cur_attr->next = malloc(sizeof(struct t_attr));
	    cur_attr = cur_attr->next;
	  }
	  cur_attr->value = NULL;
	  cur_attr->next  = NULL;
	  cur_attr->key   = calloc((s+1), sizeof(char));
	  strncpy(cur_attr->key, buf2, s);
	  s = 0;
	}
	break;
	case 6:
#line 103 "xml_scan.rl"
	{
	  struct t_attr *cur_attr = tnode[elevel]->attr;
	  while (cur_attr->next != NULL) cur_attr = cur_attr->next;
	  cur_attr->value = calloc((s+1), sizeof(char));
	  strncpy(cur_attr->value, buf2, s);
	  s = 0;
	}
	break;
	case 7:
#line 110 "xml_scan.rl"
	{
	  tnode[elevel]->body = calloc((s+1), sizeof(char));
	  strncpy(tnode[elevel]->body, buf2, s);
	  s = 0;
	}
	break;
#line 382 "xml_scan.h"
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

#line 167 "xml_scan.rl"

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
 * Auxiliary routines:
 * -------------------------------------------------------------- */

/*
 * Removes trailing and beginning whitespace together with internal
 * whitespace gaps when printing a string.
 */
void trim_printf2(size_t len, const char *str)
{
  if(len == 0) return;

  const char *end;
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;
  end++;

  int  last_ws = 1;
  char c;
  while (end > str) {
    if ((!isspace(*str)) || (last_ws == 0)) {
      c = isspace(*str) ? ' ' : *str;
      printf("%c", c);
    }
    last_ws = isspace(*str) ? 1 : 0;
    str++;
  }
}


/* --------------------------------------------------------------
 * Generic routines for accessing parsed XML tree
 * -------------------------------------------------------------- */

void print_xml_tree(const struct t_element *node, int indent) {
  struct t_element *tnode = NULL;
  struct t_attr    *tattr = NULL;
  
  if (node->body != NULL) {
    printf("%*selement %d/%s (", indent, " ", node->level, node->name);
    trim_printf2(strlen(node->body), node->body);
    printf(")\n");
  }
  else
    printf("%*selement %d/%s\n", indent, " ", node->level, node->name);

  tattr = node->attr;
  if (tattr != NULL)  printf("%*sattributes:\n", indent, " ");
  while (tattr != NULL) {
    printf("%*s  %s: ", indent, " ", tattr->key);
    trim_printf2(strlen(tattr->value), tattr->value);
    printf("\n");
    tattr = tattr->next;
  }

  tnode = node->child;
  if (tnode != NULL)  
    {
      printf("%*schildren (%s):\n", indent, " ", node->name);
      print_xml_tree(tnode,indent+2);
    }

  tnode = node->next;
  if (tnode != NULL)  
    {
      printf("%*ssiblings (%s):\n", indent, " ", node->name);
      print_xml_tree(tnode,indent);
    }
}


/*
 * Find first XML tag in tree which matches a given name and
 * (optional) different criteria checked by a given comparison
 * function.
 */
const struct t_element* find_first_element(const struct t_element *node, 
					   const  char      *name,
					   int (*cmp_function)(const struct t_element*, void*),
					   void             *context) {

  if (node == NULL) return NULL;
  struct t_element *tnode = NULL;
  const struct t_element *result = NULL;

  if (strcmp(node->name, name) == 0) {
    if (cmp_function != NULL) {
      if (cmp_function(node, context)) return node;
    }
    else 
      return node;
  }

  /* traverse children */
  tnode = node->child;
  if (tnode != NULL) {
    result = find_first_element(tnode, name, cmp_function,context);
    if (result != NULL) return result;
  }

  /* traverse siblings */
  tnode = node->next;
  if (tnode != NULL) {
    result = find_first_element(tnode, name, cmp_function,context);
    if (result != NULL) return result;
  }

  return NULL;
}


/*
 * @return first XML tag attribute matching a given name.
 */
struct t_attr* find_first_attribute(const struct t_element *node, const char* key) {
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


/*
 * @return 1, if the given grid matches the attributes
 * numberOfGridUsed, centre, subcentre, and type.
 */
int match_grid(const struct t_element *node, void *context) {
  struct t_grid_attr *ta = (struct t_grid_attr *) context;
  struct t_attr* attr = NULL;
  attr = find_first_attribute(node, "number_of_grid_used");
  if (attr == NULL) return 0;
  if (atoi(attr->value) != ta->number_of_grid_used) return 0;
  attr = find_first_attribute(node, "centre");
  if (attr == NULL) return 0;
  if (atoi(attr->value) != ta->centre) return 0;
  attr = find_first_attribute(node, "subcentre");
  if (attr == NULL) return 0;
  if (atoi(attr->value) != ta->subcentre) return 0;
  attr = find_first_attribute(node, "type");
  if (attr == NULL) return 0;
  if (strcmp(attr->value, ta->type) != 0) return 0;
  return 1;
}


/*
 * @return 1, if one of this gridset's grids matches the attributes
 * numberOfGridUsed, centre, subcentre, and type.
 */
int match_gridset(const struct t_element *node, void *context) {
  struct t_grid_attr *ta = (struct t_grid_attr *) context;
  struct t_attr*    attr  = NULL;
  struct t_element *tnode = NULL;

  /* cycle over children of type "grid" */
  tnode = node->child;
  while (tnode != NULL) {

    attr = find_first_attribute(tnode, "number_of_grid_used");
    if (attr == NULL) { tnode = tnode->next; continue; }
    if (atoi(attr->value) != ta->number_of_grid_used) { tnode = tnode->next; continue; }
    attr = find_first_attribute(tnode, "centre");
    if (attr == NULL) { tnode = tnode->next; continue; }
    if (atoi(attr->value) != ta->centre) { tnode = tnode->next; continue; }
    attr = find_first_attribute(tnode, "subcentre");
    if (attr == NULL) { tnode = tnode->next; continue; }
    if (atoi(attr->value) != ta->subcentre) { tnode = tnode->next; continue; }
    attr = find_first_attribute(tnode, "type");
    if (attr == NULL) { tnode = tnode->next; continue; }
    if (strcmp(attr->value, ta->type) != 0) { tnode = tnode->next; continue; }
    return 1;
    tnode = tnode->next;
  }
  return 0;
}


const char* get_grid_uri(struct t_grid_attr *grid_attr) {
  /* step 1: find <grid> tag: */
  const struct t_element *data_node = find_first_element(xmltree,   "data",      NULL, NULL);
  const struct t_element *list_node = find_first_element(data_node, "grid_list", NULL, NULL);
  const struct t_element *grid_node = find_first_element(list_node, "grid",      &match_grid, grid_attr);
  if (grid_node == NULL) return NULL;
  /* step 2: get corresponding URI tag: */
  const struct t_element *uri_node = find_first_element(grid_node, "uri", NULL, NULL);
  assert(uri_node       != NULL);
  assert(uri_node->body != NULL);
  return uri_node->body;
}


const char* get_grid_extpar(struct t_grid_attr *grid_attr) {
  /* step 1: find <grid> tag: */
  const struct t_element *data_node = find_first_element(xmltree,   "data",      NULL, NULL);
  const struct t_element *list_node = find_first_element(data_node, "grid_list", NULL, NULL);
  const struct t_element *grid_node = find_first_element(list_node, "grid",      &match_grid, grid_attr);
  if (grid_node == NULL) return NULL;
  /* step 2: get corresponding URI tag: */
  const struct t_element *extpar_node = find_first_element(grid_node, "extpar", NULL, NULL);
  assert(extpar_node       != NULL);
  assert(extpar_node->body != NULL);
  return extpar_node->body;
}


const char* get_grid_description(struct t_grid_attr *grid_attr) {
  /* step 1: find <grid> tag: */
  const struct t_element *data_node = find_first_element(xmltree,   "data",      NULL, NULL);
  const struct t_element *list_node = find_first_element(data_node, "grid_list", NULL, NULL);
  const struct t_element *grid_node = find_first_element(list_node, "grid",      &match_grid, grid_attr);
  if (grid_node == NULL) return NULL;
  /* step 2: get corresponding DESCR tag: */
  const struct t_element *descr_node = find_first_element(grid_node, "description", NULL, NULL);
  assert(descr_node       != NULL);
  return descr_node->body;
}

const struct t_element* find_gridset(const struct t_element* node, struct t_grid_attr *grid_attr) {
  if (node == NULL) return NULL;
  const struct t_element *gset_node = find_first_element(node, "gridset", &match_gridset, grid_attr);
  return gset_node;
}

void print_gridset(const struct t_element* gridset_node) {
  const struct t_element *descr_node = find_first_element(gridset_node, "description", NULL, NULL);
  assert(descr_node       != NULL);
  printf("- \"");
  trim_printf2(strlen(descr_node->body), descr_node->body);  
  printf("\"\n");
  
  const struct t_element *grid_node = find_first_element(gridset_node, "grid", NULL, NULL);
  while (grid_node != NULL) 
    {
      struct t_attr* attr  = NULL;
      assert(grid_node       != NULL);
      attr = find_first_attribute(grid_node, "number_of_grid_used");
      int n1 = atoi(attr->value);
      attr = find_first_attribute(grid_node, "centre");
      int n2 = atoi(attr->value);
      attr = find_first_attribute(grid_node, "subcentre");
      int n3 = atoi(attr->value);
      attr = find_first_attribute(grid_node, "type");
      printf("  numberOfGridUsed: %d, centre: %d, subcentre: %d, type: %s\n", n1,n2,n3, attr->value);  
      
      grid_node = find_first_element(grid_node->next, "grid", NULL, NULL);
    }
}

#endif // XML_SCAN_H
