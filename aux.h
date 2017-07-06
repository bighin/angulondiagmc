#ifndef __AUX_H__
#define __AUX_H__

struct vlist_t
{
	void *mem;

	size_t elementsize;
	int nelements,nalloced;
};

struct vlist_t *init_vlist(size_t elementsize,int maxsz);
void fini_vlist(struct vlist_t *lst);

void *vlist_get_element(struct vlist_t *lst,int n);
void vlist_remove_element(struct vlist_t *lst,int position);
void *vlist_add_element(struct vlist_t *lst,void *element,int position);
void *vlist_add_empty(struct vlist_t *lst,int position);
void *vlist_append(struct vlist_t *lst,void *element);
void *vlist_append_empty(struct vlist_t *lst);
int vlist_get_nr_elements(struct vlist_t *lst);

struct randomized_list_t
{
	int *items;
	int nitems,nalloced,blocksize;
};

struct randomized_list_t *init_rlist(void);
void fini_rlist(struct randomized_list_t *lst);

void rlist_add_item(struct randomized_list_t *lst,int item);
int rlist_get_random_item(struct randomized_list_t *lst);
int rlist_get_elements(struct randomized_list_t *lst);
void rlist_remove_element(struct randomized_list_t *lst,int position);
int rlist_pop_random_item(struct randomized_list_t *lst);


#endif //__AUX_H__
