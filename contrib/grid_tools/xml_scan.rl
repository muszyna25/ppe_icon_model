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

%%{
        machine xmlparse;

	action collect {
	  buf2[s] = *p; s++;
	}
	action elementStart {
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
	action elementEndSingle { tnode[elevel+1] = NULL; elevel--; s = 0; fgoto elementBody; }
	action elementEnd       { tnode[elevel+1] = NULL; elevel--; s = 0; fgoto elementBody; }
	action element          { fgoto elementBody; }
	action attributeName {
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
	action attribute {
	  struct t_attr *cur_attr = tnode[elevel]->attr;
	  while (cur_attr->next != NULL) cur_attr = cur_attr->next;
	  cur_attr->value = calloc((s+1), sizeof(char));
	  strncpy(cur_attr->value, buf2, s);
	  s = 0;
	}
	action text {
	  tnode[elevel]->body = calloc((s+1), sizeof(char));
	  strncpy(tnode[elevel]->body, buf2, s);
	  s = 0;
	}

	attribute = ^(space | [/>=])+ @collect %attributeName space* '=' space*
	  (('\'' ^'\''* @collect %attribute '\'') | ('"' ^'"'* @collect %attribute '"'));
	element = '<' space* ^(space | [/>])+ @collect %elementStart (space+ attribute)*
	  :>> (space* ((('/' %elementEndSingle)|('?' %elementEndSingle))? space* '>')  @element);
        elementBody := space* <: ((^'<'+ @collect %text) <: space*)?
	  element? :>> ('<' space* '/' ^'>'+ '>' @elementEnd);
        main := space* element space*;
}%%

%% write data nofinal;

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

  %% write init;

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
    
    %% write exec;

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
