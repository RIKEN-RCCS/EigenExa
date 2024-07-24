#pragma once

extern void eigen_libs_eigen_init_();
extern void eigen_libs_eigen_free_();
extern int eigen_blacs_eigen_get_blacs_context_();

extern void eigen_libs_eigen_sx_(const int*,const int*,double*,int*,double*,double*,int*,int*,int*,char*);
extern void eigen_libs_eigen_s_(int*,int*,double*,int*,double*,double*,int*,int*,int*,char*);

extern void eigen_libs0_eigen_get_version_(int*,char*,char*);
extern void eigen_libs0_eigen_show_version_();

extern void eigen_libs_eigen_get_matdims_(int*,int*,int*,int*,int*,char*);
extern int eigen_libs0_eigen_memory_internal_(int*,int*,int*,int*,int*);

extern void eigen_libs0_eigen_get_comm_(int*,int*,int*);
extern void eigen_libs0_eigen_get_procs_(int*,int*,int*);
extern void eigen_libs0_eigen_get_id_(int*,int*,int*);

extern int eigen_libs0_eigen_loop_start_(int*,int*,int*);
extern int eigen_libs0_eigen_loop_end_(int*,int*,int*);
extern void eigen_libs0_eigen_loop_info_(int*,int*,int*,int*,int*,int*);

extern int eigen_libs0_eigen_translate_l2g_(int*,int*,int*);
extern int eigen_libs0_eigen_translate_g2l_(int*,int*,int*);

extern int eigen_libs0_eigen_owner_node_(int*,int*,int*);
extern int eigen_libs0_eigen_owner_index_(int*,int*,int*);

extern int eigen_libs0_eigen_convert_id_xy2w_(int*,int*);
extern void eigen_libs0_eigen_convert_id_w2xy_(int*,int*,int*);

extern void eigen_libs0_eigen_get_errinfo_(int*);

typedef struct {
  double r;
  double i;
} complex;

extern void eigen_libs_eigen_h_(int*,int*,complex*,int*,double*,complex*,int*,int*,int*,char*);

