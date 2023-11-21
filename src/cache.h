#ifndef CACH_H
#define CACH_H

#define MAX_GRIDS_NUMBER 20

struct grid {
  double *x , *s_upper , *p_a, *s_upper_lower;
  double areas[3];
  int steps_number;
  double sampling_probabilities[2];
  double unif_scaler;
  double lt_properties[5];
  double rt_properties[6];
  int exist;
};

struct grids {
  struct grid grid[MAX_GRIDS_NUMBER];
  int incache;
};

extern struct grids grids;

#endif
