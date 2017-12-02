#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "rbtree_augmented.h"
#include "interval_tree_generic.h"
#include "interval_tree_double.h"

#define START(node) ((node)->start)
#define LAST(node)  ((node)->last)

INTERVAL_TREE_DEFINE(struct interval_tree_node, rb,
		     double, __subtree_last,
		     START, LAST,, interval_tree)

int interval_tree_count_overlaps(struct rb_root_cached *root, double start, double last)
{
	struct interval_tree_node *node;
	int results = 0;

	for (node = interval_tree_iter_first(root, start, last); node; node = interval_tree_iter_next(node, start, last))
		results++;

	return results;
}
