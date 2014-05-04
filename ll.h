/*
 * ll.h
 *
 * Implements common functions for maintenance of doubly linked
 * lists.  Uses macros for type a type agnostic interface.  To
 * comply with the expected structure of a doubly linked list,
 * the structure to be linked must have both a ->prev pointer
 * and a ->next pointer.
 *
 * This implementation is very simple and is intended only for
 * very simple data structures.  Do NOT attempt to put a single
 * element into two separate doubly linked list using this data
 * structure, as the ->prev and ->next pointers will be clobbered
 * between lists.
 *
 * Requires GCC extensions (typeof) for macro safety.
 */

#ifndef LL_H
#define LL_H

/*
 * Declare a doubly linked list with nodes of type _t. Usage:
 * 	LL(int) mylist = ll_new();
 * declares and initializes a new doubly linked list named my_list.
 */
#define LL(_t) struct { _t *head; _t *tail; int length; }

#define ll_head(_ll)     (_ll).head
#define ll_tail(_ll)     (_ll).tail
#define ll_length(_ll)   (_ll).length
#define ll_next(_x)      (_x)->next
#define ll_prev(_x)      (_x)->prev

#define ll_new() { NULL, NULL, 0 }

/* _n must have NULL ->prev and ->next pointers prior to calling. */
#define ll_append(_ll, _n)                     \
do {                                           \
	typeof(_n) _x = (_n);                  \
	if ((_ll).head == NULL)                \
		(_ll).head = (_ll).tail = _x;  \
	else {                                 \
		(_ll).tail->next = _x;         \
		_x->prev = (_ll).tail;         \
		(_ll).tail = _x;               \
	}                                      \
	++ll_length(_ll);                      \
} while (0)

/* _n must be in the list.  Results are undefined if it is not. */
/* (_n)'s references to the linked list are DELETED to aid in safe use. */
#define ll_delete(_ll, _n)                       \
do {                                             \
	typeof(_n) _x = (_n);                    \
	if ((_ll).head == (_ll).tail) {          \
		(_ll).head = (_ll).tail = NULL;  \
	} else if ((_ll).head == _x) {           \
		(_ll).head = (_ll).head->next;   \
		(_ll).head->prev = NULL;         \
	} else if ((_ll).tail == _x) {           \
		(_ll).tail = (_ll).tail->prev;   \
		(_ll).tail->next = NULL;         \
	} else {                                 \
		_x->prev->next = _x->next;       \
		_x->next->prev = _x->prev;       \
	}                                        \
	_x->next = _x->prev = NULL;              \
	--ll_length(_ll);                        \
} while(0)

/* Deallocate every node in the linked list using free() */
#define ll_dealloc(_ll) ll_dealloc_with((_ll), free)

/*
 * Deallocate every node in the linked list using the supplied
 * function _f.  ll_delete() ensures that _h and _t are set to
 * NULL when the last node is removed.
 */
#define ll_dealloc_with(_ll, _f)       \
do {                                   \
	typeof((_ll).head) _xx;        \
	while ((_ll).head != NULL) {   \
		_xx = (_ll).head;      \
		ll_delete((_ll), _xx); \
		(_f)(_xx);             \
	}                              \
} while (0)

#endif /* LL_H */
