#ifdef _KERNEL_CL_
#define size_t int
#define __kernel
#define __global
extern int get_global_id(int);
#endif
