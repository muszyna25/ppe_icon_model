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
	action elementEndSingle { elevel--; fgoto elementBody; }
	action elementEnd       { elevel--; fgoto elementBody; }
	action element          { fgoto elementBody; }
	action attributeName {
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
	action attribute {
	  struct t_attr *cur_attr = tnode->attr;
	  while (cur_attr->next != NULL) cur_attr = cur_attr->next;
	  cur_attr->value = malloc(sizeof(char)*(s+1));
	  strncpy(cur_attr->value, buf2, s);
	  s = 0;
	}
	action text {
	  tnode->body = malloc(sizeof(char)*(s+1));
	  strncpy(tnode->body, buf2, s);
	  s = 0;
	}

	attribute = ^(space | [/>=])+ @collect %attributeName space* '=' space*
	  (('\'' ^'\''* @collect %attribute '\'') | ('"' ^'"'* @collect %attribute '"'));
	element = '<' space* ^(space | [/>])+ @collect %elementStart (space+ attribute)*
	  :>> (space* ('/' %elementEndSingle)? space* '>' @element);
        elementBody := space* <: ((^'<'+ @collect %text) <: space*)?
	  element? :>> ('<' space* '/' ^'>'+ '>' @elementEnd);
        main := space* element space*;
}%%

%% write data nofinal;

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

#endif // XML_SCAN_H
