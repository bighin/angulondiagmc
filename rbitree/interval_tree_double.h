/* SPDX-License-Identifier: GPL-2.0 */
#ifndef _LINUX_INTERVAL_TREE_H
#define _LINUX_INTERVAL_TREE_H

#include "rbtree.h"

struct interval_tree_node {
	struct rb_node rb;
	double start;	/* Start of interval */
	double last;	/* Last location _in_ interval */
	double __subtree_last;
};

void
interval_tree_insert(struct interval_tree_node *node,
		     struct rb_root_cached *root);

void
interval_tree_remove(struct interval_tree_node *node,
		     struct rb_root_cached *root);

struct interval_tree_node *
interval_tree_iter_first(struct rb_root_cached *root,
			 double start, double last);

struct interval_tree_node *
interval_tree_iter_next(struct interval_tree_node *node,
			double start, double last);

int interval_tree_count_overlaps(struct rb_root_cached *root, double start, double last);

#endif	/* _LINUX_INTERVAL_TREE_H */
