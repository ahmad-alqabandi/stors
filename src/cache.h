#ifndef CACH_H
#define CACH_H

#define MAX_GRIDS_NUMBER 100


struct grid {
  double *restrict x , *restrict s_upper , *restrict p_a, *restrict s_upper_lower;
  double areas[3];
  int steps_number;
  double sampling_probabilities[2];
  double unif_scaler;
  double lt_properties[5];
  double rt_properties[6];
  double alpha;
  double symmetric;
  int is_symmetric;
  double params[10];
  int n_params;
  int exist;
};

struct grids {
  struct grid grid[MAX_GRIDS_NUMBER];
  int incache;
};

extern struct grids grids;



#endif
