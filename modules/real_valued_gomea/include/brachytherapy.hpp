#pragma once

#include "tools.hpp"

double dwellTimeUpperBound( int dimension );
double dwellTimeLowerBound( int dimension );

void ezilaitiniDoseDistribution();

void readReferenceFileForBrachytherapyLinearVolumeBasedInversePlanning( char *filename , int *number_of_parameters );
void ezilaitiniSetOfDoseCalculationPoints();

double getDTMRTotalViolation( double *parameters, short use_dynamic_evaluation, double scale_factor );
int* getRandomSetOfDoseCalculationPointsFromAnOrgan( int organ_id, int number_of_points );
void loadRandomSetOfDoseCalculationPoints( int *number_of_dose_calc_points );
void expandRandomSetOfDoseCalculationPoints( int *number_of_dose_calc_points );
double getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( int organ_id, int number_of_dose_calc_points_in_this_organ, int *list_of_dose_calc_points_in_this_organ, double *dose_at_point, double maximum_dose, short use_smoothening, double scale_factor );
double getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose( int organ_id, int number_of_dose_calc_points_in_this_organ, int *list_of_dose_calc_points_in_this_organ, double *dose_at_point, double minimum_dose, short use_smoothening, double scale_factor );
double getMaximumDoseInAnOrganTotalViolation( int organ_id, int number_of_dose_calc_points_in_this_organ, int *list_of_dose_calc_points_in_this_organ, double *dose_at_point, double maximum_dose );
double getMaximumDoseInAnOrgan( int organ_id, int number_of_dose_calc_points_in_this_organ, int *list_of_dose_calc_points_in_this_organ, double *dose_at_point );
void fullEvaluateBrachytherapyTreatmentPlan( double *parameters, double *objective_value, double *constraint_value );
void evaluateBrachytherapyTreatmentPlan( double *parameters, int number_of_generations, double *objective_value, double *constraint_value, short use_smoothening_volume, short use_dynamic_constraint_violation_evaluation, short use_dynamic_dtmr_constraint_evaluation );

void BT_TumorCoverage_OrganSparing_ProblemEvaluation( double *parameters, int number_of_generations, double *objective_values_result, double *constraint_values_result, short use_smoothening_volume );
void full_BT_TumorCoverage_OrganSparing_ProblemEvaluation( double *parameters, double *objective_values_result, double *constraint_values_result );
void BT_TumorCoverage_LeastOrganSparing_ProblemEvaluation( double *parameters, int number_of_generations, double *objective_values_result, double *constraint_values_result, short use_smoothening_volume );
void full_BT_TumorCoverage_LeastOrganSparing_ProblemEvaluation( double *parameters, double *objective_values_result, double *constraint_values_result );
void BT_TumorCoverage_MaximizeMinimumSafeOrgan_ProblemEvaluation( double *parameters, int number_of_generations, double *objective_values_result, double *constraint_values_result, short use_smoothening_volume );
void full_BT_TumorCoverage_MaximizeMinimumSafeOrgan_ProblemEvaluation( double *parameters, double *objective_values_result, double *constraint_values_result );
void full_BT_TumorCoverage_MaximizeMinimumSafeOrgan_ProblemEvaluation_with_Multiplier( double *parameters, double *objective_values_result, double *constraint_values_result );
void BT_TumorCoverage_MaximizeMinimumSafeOrgan_ProblemEvaluation_with_Multiplier( double *parameters, int number_of_generations, double *objective_values_result, double *constraint_values_result, short use_smoothening_volume );


void full_BT_TumorCoverage_LeastSafeIndex_ProblemEvaluation( double *parameters, double *objective_values_result, double *constraint_values_result,
                                                             double lower_bound_of_coverage, double relaxation_factor_of_safety_indices );
void BT_TumorCoverage_LeastSafeIndex_ProblemEvaluation( double *dose_distribution_buffer, double *parameters, int number_of_generations,
                                                        double *objective_values_result, double *constraint_values_result,
                                                        short use_smoothening_volume,
                                                        double lower_bound_of_coverage, double relaxation_factor_of_safety_indices );

void BT_TumorCoverage_LeastSafeIndex_PartialProblemEvaluation( double *dose_distribution_buffer, double *parameters,
                                                               double *parameters_before, int number_of_changed_dwell_times,
                                                               int *changed_dwell_time_indices, int number_of_generations,
                                                               double *objective_values_result, double *constraint_values_result,
                                                               short use_smoothening_volume,
                                                               double lower_bound_of_coverage, double relaxation_factor_of_safety_indices );

void readReferenceFileForDVHBasedInversePlanning( char *filename , char *distance_filename, int *number_of_parameters );
void ezilaitiniBrachyTherapyDVH();
double* readSamplePlan(int number_of_variables );
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#define SMOOTHENING_WITHOUT_SCALE_FACTOR 1
#define SMOOTHENING_WITH_SCALE_FACTOR 2

//#define LOWER_BOUND_OF_TREATMENT_VOLUME 0.9
//#define RELAXATION_FACTOR 2.0

double lower_bound_of_coverage;
double relaxation_factor_of_safety_indices;

int treatment_volume_id; // The ID of the organ that is treated, i.e., maximizing V100

double *lower_bound_for_dose_at_organ,
       *upper_bound_for_dose_at_organ,
       *penalty_lower_bound_at_organ,
       *penalty_upper_bound_at_organ,
       *scaled_penalty_lower_bound_at_organ,
       *scaled_penalty_upper_bound_at_organ,
       **dose_rate_matrix = NULL,
       *volume_criteria;
int    *number_of_points_in_organ,
       number_of_dwell_positions,
       number_of_dose_calculation_points,
       number_of_organ_types,
       number_of_dose_objectives,
       number_of_catheters,
       *number_of_dwell_positions_in_catheter,
       **list_of_dwell_positions_in_catheter,
       *catheter_of_dwell_position,
       *is_treatment_volume;

short  evaluated_against_full_set_of_dose_calc_points;
long   write_reference_file_now;
int    *full_number_of_points_in_organ,
       full_number_of_dose_calculation_points;
double **full_dose_rate_matrix;

int   tracking_elitist_solution;

short dwell_time_modulation_restriction,
      check_evaluation_correctness,
      dynamic_dtmr_constraint_evaluation,
      dynamic_constraint_violation_evaluation,
      use_dmax_constraints;

char *inputFilename, *distanceFilename;
short use_rounding,
      revert_to_non_rounding,
      smoothening_volume;

double brachy_tolerance = 1e-8;

int number_of_dose_distributions;
int current_individual_index;
int number_of_evaluations_per_individual;

int **list_of_dose_calc_points_in_organ,
    **full_list_of_dose_calc_points_in_organ;
double *percentage_of_dose_point_in_organ;
int **full_list_of_random_order_dose_calc_points_in_organ = NULL;

int set_of_dose_calc_points_size;
short expand_random_set_of_dose_calc_points,
      load_new_random_set_of_dose_calc_points,
      fix_random_set_of_dose_calc_points,
      load_new_data = 0;
int number_of_dvh_indices,
    *organ_id_of_dvh_index,
    *is_treatment_index,
    *is_dose_index;
double *dvh_index_dose_level,
       *dvh_index_upper_bound,
       *organ_volume;
double prescribed_dose;

double **distance_matrix;

double* readSamplePlan(int number_of_variables )
{
  double *plan, temp;
  FILE *fp;
  int i;
  plan = (double*)Malloc(number_of_variables*sizeof(double));
  fp = fopen("test_data.txt","r");
  for( i = 0; i < number_of_variables; i++)
  {
    fscanf(fp, "%lf", &temp);
    plan[i] = temp;
  }
  fclose(fp);
  return plan;
}

void BT_LeastCoverageIndex_LeastSafeIndex_PartialProblemEvaluation( double *dose_distribution_buffer, double *parameters, double *parameters_before, int number_of_changed_dwell_times, int *changed_dwell_time_indices, double *objective_values_result, double *constraint_values_result, double *dv_indices )
{
  int number_of_points_equivalent, *order, index;
  double *dose_of_points_in_this_organ, requested_volume;
  double *dose_at_point, *delta_vector, P, V, B1, B2, R1, R2, U, least_safety_index, least_coverage_index;
  double upper_bound_constraint_violation, lower_bound_constraint_violation, total_violation, v100_prostate, v150_prostate, v200_prostate, v80_vesicles, d1cc_bladder, d2cc_bladder, d1cc_rectum, d2cc_rectum, d01cc_urethra;
  int i, organ_id, point_id;

  delta_vector = (double*) Malloc(number_of_dwell_positions * sizeof( double ) );
  for( i = 0; i < number_of_dwell_positions; i++ )
    delta_vector[i] = 0;
  for( i = 0; i < number_of_changed_dwell_times; i++ )
  {
    index = changed_dwell_time_indices[i];
    delta_vector[index] = parameters[index] - parameters_before[i];
  }

  dose_at_point =  matrixVectorPartialMultiplication( dose_rate_matrix, delta_vector,
                                    number_of_dose_calculation_points, number_of_dwell_positions,
                                    number_of_changed_dwell_times, changed_dwell_time_indices );

  for( i = 0; i < number_of_dose_calculation_points; i++ )
  {
    dose_distribution_buffer[i] = dose_distribution_buffer[i] + dose_at_point[i];
    if( dose_distribution_buffer[i] < 0 )
      dose_distribution_buffer[i] = 0;
  }

  free( delta_vector );
  free( dose_at_point );
  dose_at_point = dose_distribution_buffer;

  organ_id = 1; // Prostate
  v100_prostate = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, prescribed_dose, 0, 0 );
  v150_prostate = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, 1.5*prescribed_dose, 0, 0 );
  v200_prostate = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, 2.0*prescribed_dose, 0, 0 );

  organ_id = 4; //Vesicles
  v80_vesicles = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, 0.8*prescribed_dose, 0, 0 );

  organ_id = 0; // Bladder
  dose_of_points_in_this_organ = (double*) Malloc ( number_of_points_in_organ[organ_id] * sizeof(double));
  for( i = 0; i < number_of_points_in_organ[organ_id]; i++ )
  {
    point_id = list_of_dose_calc_points_in_organ[organ_id][i];
    dose_of_points_in_this_organ[i] = dose_at_point[point_id];
  }

  order = mergeSort( dose_of_points_in_this_organ, number_of_points_in_organ[organ_id] );

  requested_volume = 1.0;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d1cc_bladder = dose_of_points_in_this_organ[ order[number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d1cc_bladder /= prescribed_dose;
  requested_volume = 2.0;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d2cc_bladder = dose_of_points_in_this_organ[ order[number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d2cc_bladder /= prescribed_dose;
  free(order);
  free(dose_of_points_in_this_organ);

  organ_id = 2; // Rectum
  dose_of_points_in_this_organ = (double*) Malloc ( number_of_points_in_organ[organ_id] * sizeof(double));
  for( i = 0; i < number_of_points_in_organ[organ_id]; i++ )
  {
    point_id = list_of_dose_calc_points_in_organ[organ_id][i];
    dose_of_points_in_this_organ[i] = dose_at_point[point_id];
  }

  order = mergeSort( dose_of_points_in_this_organ, number_of_points_in_organ[organ_id] );

  requested_volume = 1.0;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d1cc_rectum = dose_of_points_in_this_organ[ order[number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d1cc_rectum /= prescribed_dose;
  requested_volume = 2.0;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d2cc_rectum = dose_of_points_in_this_organ[ order[number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d2cc_rectum /= prescribed_dose;
  free(order);
  free(dose_of_points_in_this_organ);

  organ_id = 3; // Urethra
  dose_of_points_in_this_organ = (double*) Malloc ( number_of_points_in_organ[organ_id] * sizeof(double));
  for( i = 0; i < number_of_points_in_organ[organ_id]; i++ )
  {
    point_id = list_of_dose_calc_points_in_organ[organ_id][i];
    dose_of_points_in_this_organ[i] = dose_at_point[point_id];
  }

  order = mergeSort( dose_of_points_in_this_organ, number_of_points_in_organ[organ_id] );

  requested_volume = 0.1;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d01cc_urethra = dose_of_points_in_this_organ[ order[number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d01cc_urethra /= prescribed_dose;
  free(order);
  free(dose_of_points_in_this_organ);

  // Evaluating v150 and v200 constraints
  upper_bound_constraint_violation = 0.0;
  if( v150_prostate > 0.5 )
    upper_bound_constraint_violation += ( v150_prostate - 0.5 );
  if( v200_prostate > 0.2 )
    upper_bound_constraint_violation += ( v200_prostate - 0.2 );

  // Evaluating v100 and v80 constraints
  lower_bound_constraint_violation = 0.0;
  /*if( v100_prostate < 0.8 )
    lower_bound_constraint_violation += (0.8 - v100_prostate );
  if( v80_vesicles < 0.8 )
    lower_bound_constraint_violation += (0.8 - v80_vesicles );
  */
  // Evaluating Least Coverage Index
  P = v100_prostate - 0.95;
  V = v80_vesicles - 0.95;
  least_coverage_index = P;
  if( V < least_coverage_index )
    least_coverage_index = V;

  // Evaluating Least Safety Index
  B1 = 0.86 - d1cc_bladder;
  B2 = 0.74 - d2cc_bladder;
  R1 = 0.78 - d1cc_rectum;
  R2 = 0.74 - d2cc_rectum;
  U = 1.1 - d01cc_urethra;
  least_safety_index = B1;
  if( B2 < least_safety_index )
    least_safety_index = B2;
  if( R1 < least_safety_index )
    least_safety_index = R1;
  if( R2 < least_safety_index )
    least_safety_index = R2;
  if( U < least_safety_index )
    least_safety_index = U;

  if( least_coverage_index < -0.2 )
    lower_bound_constraint_violation += (-0.2 - least_coverage_index );
  if( least_safety_index < -0.2 )
    upper_bound_constraint_violation += (-0.2 - least_safety_index );

  *constraint_values_result = 0.0;
  total_violation = upper_bound_constraint_violation + lower_bound_constraint_violation;
  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;

  objective_values_result[0] = -least_coverage_index;
  objective_values_result[1] = -least_safety_index;

  dv_indices[0] = v100_prostate;
  dv_indices[1] = v150_prostate;
  dv_indices[2] = v200_prostate;
  dv_indices[3] = v80_vesicles;
  dv_indices[4] = d1cc_bladder;
  dv_indices[5] = d2cc_bladder;
  dv_indices[6] = d1cc_rectum;
  dv_indices[7] = d2cc_rectum;
  dv_indices[8] = d01cc_urethra;
}

void BT_LeastCoverageIndex_LeastSafeIndex_ProblemEvaluation( double *dose_distribution_buffer, double *parameters, double *objective_values_result, double *constraint_values_result, double *dv_indices )
{
  int number_of_points_equivalent, *order;
  double *dose_of_points_in_this_organ, requested_volume;
  double *dose_at_point, P, V, B1, B2, R1, R2, U, least_safety_index, least_coverage_index;
  double upper_bound_constraint_violation, lower_bound_constraint_violation, total_violation, v100_prostate, v150_prostate, v200_prostate, v80_vesicles, d1cc_bladder, d2cc_bladder, d1cc_rectum, d2cc_rectum, d01cc_urethra;
  int i, organ_id, point_id;

  dose_at_point = matrixVectorMultiplication( dose_rate_matrix, parameters, number_of_dose_calculation_points, number_of_dwell_positions );

  if( dose_distribution_buffer != NULL )
  {
    for( i = 0; i < number_of_dose_calculation_points; i++ )
      dose_distribution_buffer[i] = dose_at_point[i];
  }

  organ_id = 1; // Prostate
  v100_prostate = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, prescribed_dose, 0, 0 );
  v150_prostate = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, 1.5*prescribed_dose, 0, 0 );
  v200_prostate = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, 2.0*prescribed_dose, 0, 0 );

  organ_id = 4; //Vesicles
  v80_vesicles = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, 0.8*prescribed_dose, 0, 0 );

  organ_id = 0; // Bladder
  dose_of_points_in_this_organ = (double*) Malloc ( number_of_points_in_organ[organ_id] * sizeof(double));
  for( i = 0; i < number_of_points_in_organ[organ_id]; i++ )
  {
    point_id = list_of_dose_calc_points_in_organ[organ_id][i];
    dose_of_points_in_this_organ[i] = dose_at_point[point_id];
  }

  order = mergeSort( dose_of_points_in_this_organ, number_of_points_in_organ[organ_id] );

  requested_volume = 1.0;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d1cc_bladder = dose_of_points_in_this_organ[ order[number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d1cc_bladder /= prescribed_dose;
  requested_volume = 2.0;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d2cc_bladder = dose_of_points_in_this_organ[ order[number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d2cc_bladder /= prescribed_dose;
  free(order);
  free(dose_of_points_in_this_organ);

  organ_id = 2; // Rectum
  dose_of_points_in_this_organ = (double*) Malloc ( number_of_points_in_organ[organ_id] * sizeof(double));
  for( i = 0; i < number_of_points_in_organ[organ_id]; i++ )
  {
    point_id = list_of_dose_calc_points_in_organ[organ_id][i];
    dose_of_points_in_this_organ[i] = dose_at_point[point_id];
  }

  order = mergeSort( dose_of_points_in_this_organ, number_of_points_in_organ[organ_id] );

  requested_volume = 1.0;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d1cc_rectum = dose_of_points_in_this_organ[ order[number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d1cc_rectum /= prescribed_dose;
  requested_volume = 2.0;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d2cc_rectum = dose_of_points_in_this_organ[ order[number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d2cc_rectum /= prescribed_dose;
  free(order);
  free(dose_of_points_in_this_organ);

  organ_id = 3; // Urethra
  dose_of_points_in_this_organ = (double*) Malloc ( number_of_points_in_organ[organ_id] * sizeof(double));
  for( i = 0; i < number_of_points_in_organ[organ_id]; i++ )
  {
    point_id = list_of_dose_calc_points_in_organ[organ_id][i];
    dose_of_points_in_this_organ[i] = dose_at_point[point_id];
  }

  order = mergeSort( dose_of_points_in_this_organ, number_of_points_in_organ[organ_id] );

  requested_volume = 0.1;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d01cc_urethra = dose_of_points_in_this_organ[ order[number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d01cc_urethra /= prescribed_dose;
  free(order);
  free(dose_of_points_in_this_organ);
  free( dose_at_point );

  // Evaluating v150 and v200 constraints
  upper_bound_constraint_violation = 0.0;
  if( v150_prostate > 0.5 )
    upper_bound_constraint_violation += ( v150_prostate - 0.5 );
  if( v200_prostate > 0.2 )
    upper_bound_constraint_violation += ( v200_prostate - 0.2 );

  // Evaluating v100 and v80 constraints
  lower_bound_constraint_violation = 0.0;
  /*if( v100_prostate < 0.8 )
    lower_bound_constraint_violation += (0.8 - v100_prostate );
  if( v80_vesicles < 0.8 )
    lower_bound_constraint_violation += (0.8 - v80_vesicles );
 */
  // Evaluating Least Coverage Index
  P = v100_prostate - 0.95;
  V = v80_vesicles - 0.95;
  least_coverage_index = P;
  if( V < least_coverage_index )
    least_coverage_index = V;

  // Evaluating Least Safety Index
  B1 = 0.86 - d1cc_bladder;
  B2 = 0.74 - d2cc_bladder;
  R1 = 0.78 - d1cc_rectum;
  R2 = 0.74 - d2cc_rectum;
  U = 1.1 - d01cc_urethra;
  least_safety_index = B1;
  if( B2 < least_safety_index )
    least_safety_index = B2;
  if( R1 < least_safety_index )
    least_safety_index = R1;
  if( R2 < least_safety_index )
    least_safety_index = R2;
  if( U < least_safety_index )
    least_safety_index = U;

  if( least_coverage_index < -0.2 )
    lower_bound_constraint_violation += (-0.2 - least_coverage_index );
  if( least_safety_index < -0.2 )
    upper_bound_constraint_violation += (-0.2 - least_safety_index );

  *constraint_values_result = 0.0;
  total_violation = upper_bound_constraint_violation + lower_bound_constraint_violation;
  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;

  objective_values_result[0] = -least_coverage_index;
  objective_values_result[1] = -least_safety_index;

  dv_indices[0] = v100_prostate;
  dv_indices[1] = v150_prostate;
  dv_indices[2] = v200_prostate;
  dv_indices[3] = v80_vesicles;
  dv_indices[4] = d1cc_bladder;
  dv_indices[5] = d2cc_bladder;
  dv_indices[6] = d1cc_rectum;
  dv_indices[7] = d2cc_rectum;
  dv_indices[8] = d01cc_urethra;
}

void full_BT_LeastCoverageIndex_LeastSafeIndex_ProblemEvaluation( double *parameters, double *objective_values_result, double *constraint_values_result, double *dv_indices )
{
  int number_of_points_equivalent, *order;
  double *dose_of_points_in_this_organ, requested_volume;
  double *dose_at_point, P, V, B1, B2, R1, R2, U, least_safety_index, least_coverage_index;
  double upper_bound_constraint_violation, lower_bound_constraint_violation, total_violation, v100_prostate, v150_prostate, v200_prostate, v80_vesicles, d1cc_bladder, d2cc_bladder, d1cc_rectum, d2cc_rectum, d01cc_urethra;
  int i, organ_id, point_id;

  dose_at_point = matrixVectorMultiplication( full_dose_rate_matrix, parameters, full_number_of_dose_calculation_points, number_of_dwell_positions );

  organ_id = 1; // Prostate
  v100_prostate = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, prescribed_dose, 0, 0 );
  v150_prostate = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, 1.5*prescribed_dose, 0, 0 );
  v200_prostate = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, 2.0*prescribed_dose, 0, 0 );

  organ_id = 4; //Vesicles
  v80_vesicles = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, 0.8*prescribed_dose, 0, 0 );

  organ_id = 0; // Bladder
  dose_of_points_in_this_organ = (double*) Malloc ( full_number_of_points_in_organ[organ_id] * sizeof(double));
  for( i = 0; i < full_number_of_points_in_organ[organ_id]; i++ )
  {
    point_id = full_list_of_dose_calc_points_in_organ[organ_id][i];
    dose_of_points_in_this_organ[i] = dose_at_point[point_id];
  }

  order = mergeSort( dose_of_points_in_this_organ, full_number_of_points_in_organ[organ_id] );

  requested_volume = 1.0;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)full_number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d1cc_bladder = dose_of_points_in_this_organ[ order[full_number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d1cc_bladder /= prescribed_dose;
  requested_volume = 2.0;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)full_number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d2cc_bladder = dose_of_points_in_this_organ[ order[full_number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d2cc_bladder /= prescribed_dose;
  free(order);
  free(dose_of_points_in_this_organ);

  organ_id = 2; // Rectum
  dose_of_points_in_this_organ = (double*) Malloc ( full_number_of_points_in_organ[organ_id] * sizeof(double));
  for( i = 0; i < full_number_of_points_in_organ[organ_id]; i++ )
  {
    point_id = full_list_of_dose_calc_points_in_organ[organ_id][i];
    dose_of_points_in_this_organ[i] = dose_at_point[point_id];
  }

  order = mergeSort( dose_of_points_in_this_organ, full_number_of_points_in_organ[organ_id] );

  requested_volume = 1.0;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)full_number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d1cc_rectum = dose_of_points_in_this_organ[ order[full_number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d1cc_rectum /= prescribed_dose;
  requested_volume = 2.0;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)full_number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d2cc_rectum = dose_of_points_in_this_organ[ order[full_number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d2cc_rectum /= prescribed_dose;
  free(order);
  free(dose_of_points_in_this_organ);

  organ_id = 3; // Urethra
  dose_of_points_in_this_organ = (double*) Malloc ( full_number_of_points_in_organ[organ_id] * sizeof(double));
  for( i = 0; i < full_number_of_points_in_organ[organ_id]; i++ )
  {
    point_id = full_list_of_dose_calc_points_in_organ[organ_id][i];
    dose_of_points_in_this_organ[i] = dose_at_point[point_id];
  }

  order = mergeSort( dose_of_points_in_this_organ, full_number_of_points_in_organ[organ_id] );

  requested_volume = 0.1;
  number_of_points_equivalent = (int)( ( requested_volume * ((double)full_number_of_points_in_organ[organ_id]) ) / organ_volume[organ_id] );
  d01cc_urethra = dose_of_points_in_this_organ[ order[full_number_of_points_in_organ[organ_id] - number_of_points_equivalent] ];
  d01cc_urethra /= prescribed_dose;
  free(order);
  free(dose_of_points_in_this_organ);
  free( dose_at_point );

  // Evaluating v150 and v200 constraints
  upper_bound_constraint_violation = 0.0;
  if( v150_prostate > 0.5 )
    upper_bound_constraint_violation += ( v150_prostate - 0.5 );
  if( v200_prostate > 0.2 )
    upper_bound_constraint_violation += ( v200_prostate - 0.2 );

  // Evaluating v100 and v80 constraints
  lower_bound_constraint_violation = 0.0;
  /*if( v100_prostate < 0.8 )
    lower_bound_constraint_violation += (0.8 - v100_prostate );
  if( v80_vesicles < 0.8 )
    lower_bound_constraint_violation += (0.8 - v80_vesicles );
  */
  // Evaluating Least Coverage Index
  P = v100_prostate - 0.95;
  V = v80_vesicles - 0.95;
  least_coverage_index = P;
  if( V < least_coverage_index )
    least_coverage_index = V;

  // Evaluating Least Safety Index
  B1 = 0.86 - d1cc_bladder;
  B2 = 0.74 - d2cc_bladder;
  R1 = 0.78 - d1cc_rectum;
  R2 = 0.74 - d2cc_rectum;
  U = 1.1 - d01cc_urethra;
  least_safety_index = B1;
  if( B2 < least_safety_index )
    least_safety_index = B2;
  if( R1 < least_safety_index )
    least_safety_index = R1;
  if( R2 < least_safety_index )
    least_safety_index = R2;
  if( U < least_safety_index )
    least_safety_index = U;

  if( least_coverage_index < -0.2 )
    lower_bound_constraint_violation += (-0.2 - least_coverage_index );
  if( least_safety_index < -0.2 )
    upper_bound_constraint_violation += (-0.2 - least_safety_index );

  *constraint_values_result = 0.0;
  total_violation = upper_bound_constraint_violation + lower_bound_constraint_violation;
  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;

  objective_values_result[0] = -least_coverage_index;
  objective_values_result[1] = -least_safety_index;

  dv_indices[0] = v100_prostate;
  dv_indices[1] = v150_prostate;
  dv_indices[2] = v200_prostate;
  dv_indices[3] = v80_vesicles;
  dv_indices[4] = d1cc_bladder;
  dv_indices[5] = d2cc_bladder;
  dv_indices[6] = d1cc_rectum;
  dv_indices[7] = d2cc_rectum;
  dv_indices[8] = d01cc_urethra;
}


void BT_TumorCoverage_LeastSafeIndex_PartialProblemEvaluation( double *dose_distribution_buffer, double *parameters,
                                                               double *parameters_before, int number_of_changed_dwell_times,
                                                               int *changed_dwell_time_indices, int number_of_generations,
                                                               double *objective_values_result, double *constraint_values_result,
                                                               short use_smoothening_volume,
                                                               double lower_bound_of_coverage, double relaxation_factor_of_safety_indices )
{
  double *dose_at_point, *delta_vector;
  double tumor_coverage, volume, safe_index, least_treatment_index_value, least_safe_index_value, upper_bound_constraint_violation, lower_bound_constraint_violation, dtmr_constraint_violation, total_violation;
  int i, index, organ_id, count_tumor, count_organ, dvh_index;

  count_tumor = 0;
  count_organ = 0;
  tumor_coverage = 0.0;
  upper_bound_constraint_violation = 0.0;
  lower_bound_constraint_violation = 0.0;
  dtmr_constraint_violation = 0.0;
  least_treatment_index_value = 1.0;
  least_safe_index_value = 1.0;

  delta_vector = (double*) Malloc(number_of_dwell_positions * sizeof( double ) );
  for( i = 0; i < number_of_dwell_positions; i++ )
    delta_vector[i] = 0;
  for( i = 0; i < number_of_changed_dwell_times; i++ )
  {
    index = changed_dwell_time_indices[i];
    delta_vector[index] = parameters[index] - parameters_before[i];
  }

  dose_at_point =  matrixVectorPartialMultiplication( dose_rate_matrix, delta_vector,
                                    number_of_dose_calculation_points, number_of_dwell_positions,
                                    number_of_changed_dwell_times, changed_dwell_time_indices );

  for( i = 0; i < number_of_dose_calculation_points; i++ )
  {
    dose_distribution_buffer[i] = dose_distribution_buffer[i] + dose_at_point[i];
    if( dose_distribution_buffer[i] < 0 )
      dose_distribution_buffer[i] = 0;
  }

  free( delta_vector );
  free( dose_at_point );
  dose_at_point = dose_distribution_buffer;

  if ( dwell_time_modulation_restriction )
    dtmr_constraint_violation = getDTMRTotalViolation( parameters, dynamic_dtmr_constraint_evaluation, (double)number_of_generations );
  else
    dtmr_constraint_violation = 0.0;

  for( dvh_index = 0; dvh_index < number_of_dvh_indices; dvh_index++ )
  {
    organ_id = organ_id_of_dvh_index[dvh_index];
    if( is_treatment_index[dvh_index] )
    {
      if( use_smoothening_volume )
        volume = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, dvh_index_dose_level[dvh_index], SMOOTHENING_WITHOUT_SCALE_FACTOR, 0 );
      else
        volume = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, dvh_index_dose_level[dvh_index], 0, 0 );

      if( volume < lower_bound_of_coverage )
      {
        lower_bound_constraint_violation += ( lower_bound_of_coverage - volume );
        tumor_coverage = 0.0;
      }
      else
      {
//        tumor_coverage = ( volume - LOWER_BOUND_OF_TREATMENT_VOLUME ) / ( 1.0 - LOWER_BOUND_OF_TREATMENT_VOLUME );
        tumor_coverage = volume;
      }

      if( tumor_coverage < least_treatment_index_value )
        least_treatment_index_value = tumor_coverage;

      count_tumor++;
    }
    else
    {
      volume = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, dvh_index_dose_level[dvh_index], 0, 0 );
      if( volume > relaxation_factor_of_safety_indices*dvh_index_upper_bound[dvh_index] )
      {
        upper_bound_constraint_violation += (volume - relaxation_factor_of_safety_indices*dvh_index_upper_bound[dvh_index]);
        safe_index = 0.0;
      }
      else
      {
        safe_index = 1.0 - (volume / (relaxation_factor_of_safety_indices*dvh_index_upper_bound[dvh_index]));
      }
      //printf("organ_id: %d - dose_level: %lf - volume: %lf upper bound: %lf\n", organ_id, dvh_index_dose_level[dvh_index], volume, dvh_index_upper_bound[dvh_index]);

      if( safe_index < least_safe_index_value )
        least_safe_index_value = safe_index;
      count_organ++;
    }
  }
  if( count_tumor > 0)
    objective_values_result[0] = -least_treatment_index_value;
  else
    objective_values_result[0] = 0;
  if( count_organ >0 )
    objective_values_result[1] = -least_safe_index_value;
  else
    objective_values_result[1] = 0;

  total_violation = dtmr_constraint_violation + upper_bound_constraint_violation + lower_bound_constraint_violation;

  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;
  else
    *constraint_values_result = 0;
}

void BT_TumorCoverage_LeastSafeIndex_ProblemEvaluation( double *dose_distribution_buffer, double *parameters, int number_of_generations,
                                                        double *objective_values_result, double *constraint_values_result,
                                                        short use_smoothening_volume,
                                                        double lower_bound_of_coverage, double relaxation_factor_of_safety_indices )
{
  double *dose_at_point;
  double tumor_coverage, volume, safe_index, least_treatment_index_value, least_safe_index_value, upper_bound_constraint_violation, lower_bound_constraint_violation, dtmr_constraint_violation, total_violation;
  int i, organ_id, count_tumor, count_organ, dvh_index;

  count_tumor = 0;
  count_organ = 0;
  tumor_coverage = 0.0;
  upper_bound_constraint_violation = 0.0;
  lower_bound_constraint_violation = 0.0;
  dtmr_constraint_violation = 0.0;
  least_treatment_index_value = 1.0;
  least_safe_index_value = 1.0;

  dose_at_point = matrixVectorMultiplication( dose_rate_matrix, parameters, number_of_dose_calculation_points, number_of_dwell_positions );

  if( dose_distribution_buffer != NULL )
  {
    for( i = 0; i < number_of_dose_calculation_points; i++ )
      dose_distribution_buffer[i] = dose_at_point[i];
  }

  if ( dwell_time_modulation_restriction )
    dtmr_constraint_violation = getDTMRTotalViolation( parameters, dynamic_dtmr_constraint_evaluation, (double)number_of_generations );
  else
    dtmr_constraint_violation = 0.0;

  for( dvh_index = 0; dvh_index < number_of_dvh_indices; dvh_index++ )
  {
    organ_id = organ_id_of_dvh_index[dvh_index];
    if( is_treatment_index[dvh_index] )
    {
      if( use_smoothening_volume )
        volume = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, dvh_index_dose_level[dvh_index], SMOOTHENING_WITHOUT_SCALE_FACTOR, 0 );
      else
        volume = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, dvh_index_dose_level[dvh_index], 0, 0 );

      if( volume < lower_bound_of_coverage )
      {
        lower_bound_constraint_violation += ( lower_bound_of_coverage - volume );
        tumor_coverage = 0.0;
      }
      else
      {
//        tumor_coverage = ( volume - LOWER_BOUND_OF_TREATMENT_VOLUME ) / ( 1.0 - LOWER_BOUND_OF_TREATMENT_VOLUME );
        tumor_coverage = volume;
      }

      if( tumor_coverage < least_treatment_index_value )
        least_treatment_index_value = tumor_coverage;

      count_tumor++;
    }
    else
    {
      volume = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, dvh_index_dose_level[dvh_index], 0, 0 );
      if( volume > relaxation_factor_of_safety_indices*dvh_index_upper_bound[dvh_index] )
      {
        upper_bound_constraint_violation += (volume - relaxation_factor_of_safety_indices*dvh_index_upper_bound[dvh_index]);
        safe_index = 0.0;
      }
      else
      {
        safe_index = 1.0 - (volume / (relaxation_factor_of_safety_indices*dvh_index_upper_bound[dvh_index]));
      }
      //printf("organ_id: %d - dose_level: %lf - volume: %lf upper bound: %lf\n", organ_id, dvh_index_dose_level[dvh_index], volume, dvh_index_upper_bound[dvh_index]);

      if( safe_index < least_safe_index_value )
        least_safe_index_value = safe_index;
      count_organ++;
    }
  }
  if( count_tumor > 0)
    objective_values_result[0] = -least_treatment_index_value;
  else
    objective_values_result[0] = 0;
  if( count_organ >0 )
    objective_values_result[1] = -least_safe_index_value;
  else
    objective_values_result[1] = 0;

  total_violation = dtmr_constraint_violation + upper_bound_constraint_violation + lower_bound_constraint_violation;

  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;
  else
    *constraint_values_result = 0;

  free( dose_at_point );
}

void full_BT_TumorCoverage_LeastSafeIndex_ProblemEvaluation( double *parameters, double *objective_values_result, double *constraint_values_result,
                                                             double lower_bound_of_coverage, double relaxation_factor_of_safety_indices )
{
  double *dose_at_point;
  double tumor_coverage, volume, safe_index, least_treatment_index_value, least_safe_index_value, upper_bound_constraint_violation, lower_bound_constraint_violation, dtmr_constraint_violation, total_violation;
  int organ_id, count_tumor, count_organ, dvh_index;

  count_tumor = 0;
  count_organ = 0;
  tumor_coverage = 0.0;
  upper_bound_constraint_violation = 0.0;
  lower_bound_constraint_violation = 0.0;
  dtmr_constraint_violation = 0.0;
  least_treatment_index_value = 1.0;
  least_safe_index_value = 1.0;

  dose_at_point = matrixVectorMultiplication( full_dose_rate_matrix, parameters, full_number_of_dose_calculation_points, number_of_dwell_positions );

  if ( dwell_time_modulation_restriction )
    dtmr_constraint_violation = getDTMRTotalViolation( parameters, 0, 0 );
  else
    dtmr_constraint_violation = 0.0;

  for( dvh_index = 0; dvh_index < number_of_dvh_indices; dvh_index++ )
  {
    organ_id = organ_id_of_dvh_index[dvh_index];

    if( is_treatment_index[dvh_index] )
    {
      volume = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, dvh_index_dose_level[dvh_index], 0, 0 );

      if( volume < lower_bound_of_coverage )
      {
        lower_bound_constraint_violation += ( lower_bound_of_coverage - volume );
        tumor_coverage = 0.0;
      }
      else
      {
//        tumor_coverage = ( volume - LOWER_BOUND_OF_TREATMENT_VOLUME ) / ( 1.0 - LOWER_BOUND_OF_TREATMENT_VOLUME );
        tumor_coverage = volume;
      }

      if( tumor_coverage < least_treatment_index_value )
        least_treatment_index_value = tumor_coverage;

      count_tumor++;
    }
    else
    {
      volume = getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose( organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, dvh_index_dose_level[dvh_index], 0, 0 );
      if( volume > relaxation_factor_of_safety_indices*dvh_index_upper_bound[dvh_index] )
      {
        upper_bound_constraint_violation += (volume - relaxation_factor_of_safety_indices*dvh_index_upper_bound[dvh_index]);
        safe_index = 0.0;
      }
      else
      {
        safe_index = 1.0 - (volume / (relaxation_factor_of_safety_indices*dvh_index_upper_bound[dvh_index]));
      }
      if( safe_index < least_safe_index_value )
        least_safe_index_value = safe_index;
      count_organ++;
    }
  }
  if( count_tumor > 0)
    objective_values_result[0] = -least_treatment_index_value;
  else
    objective_values_result[0] = 0;
  if( count_organ >0 )
    objective_values_result[1] = -least_safe_index_value;
  else
    objective_values_result[1] = 0;

  total_violation = dtmr_constraint_violation + upper_bound_constraint_violation + lower_bound_constraint_violation;

  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;
  else
    *constraint_values_result = 0;

  free( dose_at_point );
}


void BT_TumorCoverage_MaximizeMinimumSafeOrgan_ProblemEvaluation_with_Multiplier( double *parameters, int number_of_generations, double *objective_values_result, double *constraint_values_result, short use_smoothening_volume )
{
  double *dose_at_point, *dwell_times;
  double tumor_coverage, organ_sparing, minimum_obj, upper_bound_constraint_violation, dtmr_constraint_violation, total_violation, dmax, converted_dmax;
  int organ_id, count_tumor, count_organ, i;

  count_tumor = 0;
  count_organ = 0;
  tumor_coverage = 0.0;
  organ_sparing = 0.0;
  minimum_obj = 1.0;
  upper_bound_constraint_violation = 0.0;
  dtmr_constraint_violation = 0.0;

  dwell_times = (double*) Malloc( number_of_dwell_positions * sizeof(double) );
  for( i = 0; i < number_of_dwell_positions; i++)
  {
    dwell_times[i] = parameters[i]*parameters[number_of_dwell_positions];
  }

  dose_at_point = matrixVectorMultiplication( dose_rate_matrix, dwell_times, number_of_dose_calculation_points, number_of_dwell_positions );

  if ( dwell_time_modulation_restriction )
    dtmr_constraint_violation = getDTMRTotalViolation( dwell_times, dynamic_dtmr_constraint_evaluation, (double)number_of_generations );
  else
    dtmr_constraint_violation = 0.0;

  for( organ_id = 0; organ_id < number_of_organ_types; organ_id++ )
  {
    if( is_treatment_volume[organ_id] )
    {
      if( use_smoothening_volume )
        tumor_coverage += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], SMOOTHENING_WITHOUT_SCALE_FACTOR, 0 );
      else
        tumor_coverage += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );
      count_tumor++;
    }
    else
    {
      dmax = getMaximumDoseInAnOrgan( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point );
      if ( dmax > 2.0*upper_bound_for_dose_at_organ[organ_id] )
      {
        //upper_bound_constraint_violation += (2.0*upper_bound_for_dose_at_organ[organ_id] - dmax);
        upper_bound_constraint_violation += getMaximumDoseInAnOrganTotalViolation( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, 2.0*upper_bound_for_dose_at_organ[organ_id]);
        converted_dmax = 0;
      }
      else
        converted_dmax = 1.0 - (dmax/(2.0*upper_bound_for_dose_at_organ[organ_id]));
      if( converted_dmax < minimum_obj )
        minimum_obj = converted_dmax;

      if( use_smoothening_volume )
        organ_sparing = getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], SMOOTHENING_WITH_SCALE_FACTOR, (double) number_of_generations );
      else
        organ_sparing = getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );
      if( organ_sparing < minimum_obj )
        minimum_obj = organ_sparing;
      count_organ++;

    }
  }
  if( count_tumor > 0)
    objective_values_result[0] = -(tumor_coverage/count_tumor);
  else
    objective_values_result[0] = 0;
  if( count_organ >0 )
    objective_values_result[1] = -minimum_obj;
  else
    objective_values_result[1] = 0;

  total_violation = dtmr_constraint_violation + upper_bound_constraint_violation;

  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;
  else
    *constraint_values_result = 0;

  free( dose_at_point );
  free( dwell_times );
}

void full_BT_TumorCoverage_MaximizeMinimumSafeOrgan_ProblemEvaluation_with_Multiplier( double *parameters, double *objective_values_result, double *constraint_values_result )
{
  double *dose_at_point, *dwell_times;
  double tumor_coverage, organ_sparing, minimum_obj, upper_bound_constraint_violation, dtmr_constraint_violation, total_violation, dmax, converted_dmax;
  int organ_id, count_tumor, count_organ, i;

  count_tumor = 0;
  count_organ = 0;
  tumor_coverage = 0.0;
  organ_sparing = 0.0;
  minimum_obj = 1.0;
  upper_bound_constraint_violation = 0.0;
  dtmr_constraint_violation = 0.0;

  dwell_times = (double*) Malloc( number_of_dwell_positions * sizeof(double) );
  for( i = 0; i < number_of_dwell_positions; i++)
    dwell_times[i] = parameters[i]*parameters[number_of_dwell_positions];


  dose_at_point = matrixVectorMultiplication( full_dose_rate_matrix, dwell_times, full_number_of_dose_calculation_points, number_of_dwell_positions );

  if ( dwell_time_modulation_restriction )
    dtmr_constraint_violation = getDTMRTotalViolation( dwell_times, 0, 0 );
  else
    dtmr_constraint_violation = 0.0;

  for( organ_id = 0; organ_id < number_of_organ_types; organ_id++ )
  {
    if( is_treatment_volume[organ_id] )
    {
      tumor_coverage += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );
      count_tumor++;
    }
    else
    {
      dmax = getMaximumDoseInAnOrgan( organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point );
      if( dmax > 2.0*upper_bound_for_dose_at_organ[organ_id] )
      {
        //upper_bound_constraint_violation += (2.0*upper_bound_for_dose_at_organ[organ_id] - dmax);
        upper_bound_constraint_violation += getMaximumDoseInAnOrganTotalViolation( organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, 2.0*upper_bound_for_dose_at_organ[organ_id]);
        converted_dmax = 0;
      }
      else
        converted_dmax = 1 - (dmax/(2.0*upper_bound_for_dose_at_organ[organ_id]));
      if( converted_dmax < minimum_obj )
        minimum_obj = converted_dmax;

      organ_sparing = getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );

      if( organ_sparing < minimum_obj)
        minimum_obj = organ_sparing;

      count_organ++;

    }
  }
  if( count_tumor > 0)
    objective_values_result[0] = -(tumor_coverage/count_tumor);
  else
    objective_values_result[0] = 0;
  if( count_organ >0 )
    objective_values_result[1] = -minimum_obj;
  else
    objective_values_result[1] = 0;

  total_violation = dtmr_constraint_violation + upper_bound_constraint_violation;

  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;
  else
    *constraint_values_result = 0;

  free( dose_at_point );
  free( dwell_times );
}


void BT_TumorCoverage_MaximizeMinimumSafeOrgan_ProblemEvaluation( double *parameters, int number_of_generations, double *objective_values_result, double *constraint_values_result, short use_smoothening_volume )
{
  double *dose_at_point;
  double tumor_coverage, organ_sparing, minimum_obj, upper_bound_constraint_violation, dtmr_constraint_violation, total_violation, dmax, converted_dmax;
  int organ_id, count_tumor, count_organ;

  count_tumor = 0;
  count_organ = 0;
  tumor_coverage = 0.0;
  organ_sparing = 0.0;
  minimum_obj = 1.0;
  upper_bound_constraint_violation = 0.0;
  dtmr_constraint_violation = 0.0;

  dose_at_point = matrixVectorMultiplication( dose_rate_matrix, parameters, number_of_dose_calculation_points, number_of_dwell_positions );

  if ( dwell_time_modulation_restriction )
    dtmr_constraint_violation = getDTMRTotalViolation( parameters, dynamic_dtmr_constraint_evaluation, (double)number_of_generations );
  else
    dtmr_constraint_violation = 0.0;

  for( organ_id = 0; organ_id < number_of_organ_types; organ_id++ )
  {
    if( is_treatment_volume[organ_id] )
    {
      if( use_smoothening_volume )
        tumor_coverage += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], SMOOTHENING_WITHOUT_SCALE_FACTOR, 0 );
      else
        tumor_coverage += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );
      count_tumor++;
    }
    else
    {
      dmax = getMaximumDoseInAnOrgan( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point );
      if ( dmax > 2.0*upper_bound_for_dose_at_organ[organ_id] )
      {
        //upper_bound_constraint_violation += (2.0*upper_bound_for_dose_at_organ[organ_id] - dmax);
        upper_bound_constraint_violation += getMaximumDoseInAnOrganTotalViolation( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, 2.0*upper_bound_for_dose_at_organ[organ_id]);
        converted_dmax = 0;
      }
      else
        converted_dmax = 1.0 - (dmax/(2.0*upper_bound_for_dose_at_organ[organ_id]));
      if( converted_dmax < minimum_obj )
        minimum_obj = converted_dmax;

      if( use_smoothening_volume )
        organ_sparing = getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], SMOOTHENING_WITH_SCALE_FACTOR, (double) number_of_generations );
      else
        organ_sparing = getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );
      if( organ_sparing < minimum_obj )
        minimum_obj = organ_sparing;
      count_organ++;

    }
  }
  if( count_tumor > 0)
    objective_values_result[0] = -(tumor_coverage/count_tumor);
  else
    objective_values_result[0] = 0;
  if( count_organ >0 )
    objective_values_result[1] = -minimum_obj;
  else
    objective_values_result[1] = 0;

  total_violation = dtmr_constraint_violation + upper_bound_constraint_violation;

  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;
  else
    *constraint_values_result = 0;

  free( dose_at_point );
}

void full_BT_TumorCoverage_MaximizeMinimumSafeOrgan_ProblemEvaluation( double *parameters, double *objective_values_result, double *constraint_values_result )
{
  double *dose_at_point;
  double tumor_coverage, organ_sparing, minimum_obj, upper_bound_constraint_violation, dtmr_constraint_violation, total_violation, dmax, converted_dmax;
  int organ_id, count_tumor, count_organ;

  count_tumor = 0;
  count_organ = 0;
  tumor_coverage = 0.0;
  organ_sparing = 0.0;
  minimum_obj = 1.0;
  upper_bound_constraint_violation = 0.0;
  dtmr_constraint_violation = 0.0;

  dose_at_point = matrixVectorMultiplication( full_dose_rate_matrix, parameters, full_number_of_dose_calculation_points, number_of_dwell_positions );

  if ( dwell_time_modulation_restriction )
    dtmr_constraint_violation = getDTMRTotalViolation( parameters, 0, 0 );
  else
    dtmr_constraint_violation = 0.0;

  for( organ_id = 0; organ_id < number_of_organ_types; organ_id++ )
  {
    if( is_treatment_volume[organ_id] )
    {
      tumor_coverage += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );
      count_tumor++;
    }
    else
    {
      dmax = getMaximumDoseInAnOrgan( organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point );
      if( dmax > 2.0*upper_bound_for_dose_at_organ[organ_id] )
      {
        //upper_bound_constraint_violation += (2.0*upper_bound_for_dose_at_organ[organ_id] - dmax);
        upper_bound_constraint_violation += getMaximumDoseInAnOrganTotalViolation( organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, 2.0*upper_bound_for_dose_at_organ[organ_id]);
        converted_dmax = 0;
      }
      else
        converted_dmax = 1 - (dmax/(2.0*upper_bound_for_dose_at_organ[organ_id]));
      if( converted_dmax < minimum_obj )
        minimum_obj = converted_dmax;

      organ_sparing = getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );

      if( organ_sparing < minimum_obj)
        minimum_obj = organ_sparing;

      count_organ++;

    }
  }
  if( count_tumor > 0)
    objective_values_result[0] = -(tumor_coverage/count_tumor);
  else
    objective_values_result[0] = 0;
  if( count_organ >0 )
    objective_values_result[1] = -minimum_obj;
  else
    objective_values_result[1] = 0;

  total_violation = dtmr_constraint_violation + upper_bound_constraint_violation;

  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;
  else
    *constraint_values_result = 0;

  free( dose_at_point );
}


void BT_TumorCoverage_LeastOrganSparing_ProblemEvaluation( double *parameters, int number_of_generations, double *objective_values_result, double *constraint_values_result, short use_smoothening_volume )
{
  double *dose_at_point;
  double tumor_coverage, organ_sparing, least_organ_sparing, upper_bound_constraint_violation, dtmr_constraint_violation, total_violation;
  int organ_id, count_tumor, count_organ;

  count_tumor = 0;
  count_organ = 0;
  tumor_coverage = 0.0;
  organ_sparing = 0.0;
  least_organ_sparing = 1.0;
  upper_bound_constraint_violation = 0.0;
  dtmr_constraint_violation = 0.0;

  dose_at_point = matrixVectorMultiplication( dose_rate_matrix, parameters, number_of_dose_calculation_points, number_of_dwell_positions );

  if ( dwell_time_modulation_restriction )
    dtmr_constraint_violation = getDTMRTotalViolation( parameters, dynamic_dtmr_constraint_evaluation, (double)number_of_generations );
  else
    dtmr_constraint_violation = 0.0;

  for( organ_id = 0; organ_id < number_of_organ_types; organ_id++ )
  {
    if( is_treatment_volume[organ_id] )
    {
      if( use_smoothening_volume )
        tumor_coverage += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], SMOOTHENING_WITHOUT_SCALE_FACTOR, 0 );
      else
        tumor_coverage += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );
      count_tumor++;
    }
    else
    {
      if( use_dmax_constraints )
        upper_bound_constraint_violation += getMaximumDoseInAnOrganTotalViolation( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, upper_bound_for_dose_at_organ[organ_id]);

      if( use_smoothening_volume )
        organ_sparing = getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], SMOOTHENING_WITH_SCALE_FACTOR, (double) number_of_generations );
      else
        organ_sparing = getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );
      if( organ_sparing < least_organ_sparing )
        least_organ_sparing = organ_sparing;
      count_organ++;

    }
  }
  if( count_tumor > 0)
    objective_values_result[0] = -(tumor_coverage/count_tumor);
  else
    objective_values_result[0] = 0;
  if( count_organ >0 )
    objective_values_result[1] = -least_organ_sparing;
  else
    objective_values_result[1] = 0;

  total_violation = dtmr_constraint_violation + upper_bound_constraint_violation;

  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;
  else
    *constraint_values_result = 0;

  free( dose_at_point );
}

void full_BT_TumorCoverage_LeastOrganSparing_ProblemEvaluation( double *parameters, double *objective_values_result, double *constraint_values_result )
{
  double *dose_at_point;
  double tumor_coverage, organ_sparing, least_organ_sparing, upper_bound_constraint_violation, dtmr_constraint_violation, total_violation;
  int organ_id, count_tumor, count_organ;

  count_tumor = 0;
  count_organ = 0;
  tumor_coverage = 0.0;
  organ_sparing = 0.0;
  least_organ_sparing = 1.0;
  upper_bound_constraint_violation = 0.0;
  dtmr_constraint_violation = 0.0;

  dose_at_point = matrixVectorMultiplication( full_dose_rate_matrix, parameters, full_number_of_dose_calculation_points, number_of_dwell_positions );

  if ( dwell_time_modulation_restriction )
    dtmr_constraint_violation = getDTMRTotalViolation( parameters, 0, 0 );
  else
    dtmr_constraint_violation = 0.0;

  for( organ_id = 0; organ_id < number_of_organ_types; organ_id++ )
  {
    if( is_treatment_volume[organ_id] )
    {
      tumor_coverage += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );
      count_tumor++;
    }
    else
    {
      if( use_dmax_constraints )
        upper_bound_constraint_violation += getMaximumDoseInAnOrganTotalViolation( organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, upper_bound_for_dose_at_organ[organ_id]);

      organ_sparing = getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );

      if( organ_sparing < least_organ_sparing )
        least_organ_sparing = organ_sparing;

      count_organ++;

    }
  }
  if( count_tumor > 0)
    objective_values_result[0] = -(tumor_coverage/count_tumor);
  else
    objective_values_result[0] = 0;
  if( count_organ >0 )
    objective_values_result[1] = -least_organ_sparing;
  else
    objective_values_result[1] = 0;

  total_violation = dtmr_constraint_violation + upper_bound_constraint_violation;

  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;
  else
    *constraint_values_result = 0;

  free( dose_at_point );
}


void BT_TumorCoverage_OrganSparing_ProblemEvaluation( double *parameters, int number_of_generations, double *objective_values_result, double *constraint_values_result, short use_smoothening_volume )
{
  double *dose_at_point;
  double tumor_coverage, organ_sparing, upper_bound_constraint_violation, dtmr_constraint_violation, total_violation;
  int organ_id, count_tumor, count_organ;

  count_tumor = 0;
  count_organ = 0;
  tumor_coverage = 0.0;
  organ_sparing = 0.0;
  upper_bound_constraint_violation = 0.0;
  dtmr_constraint_violation = 0.0;

  dose_at_point = matrixVectorMultiplication( dose_rate_matrix, parameters, number_of_dose_calculation_points, number_of_dwell_positions );

  if ( dwell_time_modulation_restriction )
    dtmr_constraint_violation = getDTMRTotalViolation( parameters, dynamic_dtmr_constraint_evaluation, (double)number_of_generations );
  else
    dtmr_constraint_violation = 0.0;

  for( organ_id = 0; organ_id < number_of_organ_types; organ_id++ )
  {
    if( is_treatment_volume[organ_id] )
    {
      if( use_smoothening_volume )
        tumor_coverage += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], SMOOTHENING_WITHOUT_SCALE_FACTOR, 0 );
      else
        tumor_coverage += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );
      count_tumor++;
    }
    else
    {
      if( use_dmax_constraints )
        upper_bound_constraint_violation += getMaximumDoseInAnOrganTotalViolation( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, upper_bound_for_dose_at_organ[organ_id]);

      if( use_smoothening_volume )
        organ_sparing += getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], SMOOTHENING_WITH_SCALE_FACTOR, (double) number_of_generations );
      else
        organ_sparing += getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );
      count_organ++;

    }
  }
  if( count_tumor > 0)
    objective_values_result[0] = -(tumor_coverage/count_tumor);
  else
    objective_values_result[0] = 0;
  if( count_organ >0 )
    objective_values_result[1] = -(organ_sparing/count_organ);
  else
    objective_values_result[1] = 0;

  total_violation = dtmr_constraint_violation + upper_bound_constraint_violation;

  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;
  else
    *constraint_values_result = 0;

  free( dose_at_point );
}

void full_BT_TumorCoverage_OrganSparing_ProblemEvaluation( double *parameters, double *objective_values_result, double *constraint_values_result )
{
  double *dose_at_point;
  double tumor_coverage, organ_sparing, upper_bound_constraint_violation, dtmr_constraint_violation, total_violation;
  int organ_id, count_tumor, count_organ;

  count_tumor = 0;
  count_organ = 0;
  tumor_coverage = 0.0;
  organ_sparing = 0.0;
  upper_bound_constraint_violation = 0.0;
  dtmr_constraint_violation = 0.0;

  dose_at_point = matrixVectorMultiplication( full_dose_rate_matrix, parameters, full_number_of_dose_calculation_points, number_of_dwell_positions );

  if ( dwell_time_modulation_restriction )
    dtmr_constraint_violation = getDTMRTotalViolation( parameters, 0, 0 );
  else
    dtmr_constraint_violation = 0.0;

  for( organ_id = 0; organ_id < number_of_organ_types; organ_id++ )
  {
    if( is_treatment_volume[organ_id] )
    {
      tumor_coverage += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );
      count_tumor++;
    }
    else
    {
      if( use_dmax_constraints )
        upper_bound_constraint_violation += getMaximumDoseInAnOrganTotalViolation( organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, upper_bound_for_dose_at_organ[organ_id]);

      organ_sparing += getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );

      count_organ++;

    }
  }
  if( count_tumor > 0)
    objective_values_result[0] = -(tumor_coverage/count_tumor);
  else
    objective_values_result[0] = 0;
  if( count_organ >0 )
    objective_values_result[1] = -(organ_sparing/count_organ);
  else
    objective_values_result[1] = 0;

  total_violation = dtmr_constraint_violation + upper_bound_constraint_violation;

  if( total_violation > brachy_tolerance )
    *constraint_values_result = total_violation;
  else
    *constraint_values_result = 0;

  free( dose_at_point );
}

double dwellTimeUpperBound( int dimension )
{
  return ( 60.0 ); // 60 seconds
}

double dwellTimeLowerBound( int dimension )
{
  return ( 0.0 );
}

void evaluateBrachytherapyTreatmentPlan( double *parameters, int number_of_generations, double *objective_value, double *constraint_value, short use_smoothening_volume, short use_dynamic_constraint_violation_evaluation, short use_dynamic_dtmr_constraint_evaluation )
{
  double *dose_at_point, total_violation, obj, volume;
  int     i, organ_id;
  double *backup = NULL;
  double  penalty_upper_bound_constraint,
          penalty_volume_constraint,
          penalty_dtmr_constraint;
  double  upper_bound_constraint_violation,
          volume_constraint_violation,
          dtmr_constraint_violation;

  if(use_rounding)
  {
    backup = (double*)Malloc(number_of_dwell_positions*sizeof(double));
    for( i = 0; i < number_of_dwell_positions; i++ )
    {
      backup[i] = parameters[i];
      parameters[i] = (long int)( parameters[i] * 10 + 0.5 ) / 10.0;
    }
  }

  dose_at_point = matrixVectorMultiplication( dose_rate_matrix, parameters, number_of_dose_calculation_points, number_of_dwell_positions );

  if ( dwell_time_modulation_restriction )
    dtmr_constraint_violation = getDTMRTotalViolation( parameters, use_dynamic_dtmr_constraint_evaluation, (double)number_of_generations );
  else
    dtmr_constraint_violation = 0.0;

  obj = 0.0;
  volume_constraint_violation = 0.0;
  upper_bound_constraint_violation = 0.0;

  for( organ_id = 0; organ_id < number_of_organ_types; organ_id++)
  {
    if( is_treatment_volume[organ_id] )
    {
      if( use_smoothening_volume )
        obj += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], SMOOTHENING_WITHOUT_SCALE_FACTOR, 0 );
      else
        obj += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0 , 0 );
    }
    else
    {
      if( use_smoothening_volume )
        volume = getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], SMOOTHENING_WITH_SCALE_FACTOR, (double) number_of_generations );
      else
        volume = getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );

      if( (volume < volume_criteria[organ_id]) && (fabs(volume - volume_criteria[organ_id] ) >  brachy_tolerance ))
        volume_constraint_violation += ( volume_criteria[organ_id] - volume);
     //printf("%d: %lf - %lf\n", organ_id, volume_criteria[organ_id], volume);
      if( use_dmax_constraints )
        upper_bound_constraint_violation += getMaximumDoseInAnOrganTotalViolation( organ_id, number_of_points_in_organ[organ_id], list_of_dose_calc_points_in_organ[organ_id], dose_at_point, upper_bound_for_dose_at_organ[organ_id]);
    }
  }

  // Convert maximization of V100_PTV to minimization
  *objective_value = -obj;

  // Process constraint violations
  if( use_dynamic_constraint_violation_evaluation )
  {
    if(number_of_generations > 0 )
    {
      penalty_upper_bound_constraint = number_of_generations;
      penalty_volume_constraint = number_of_generations;
      penalty_dtmr_constraint = number_of_generations;
    }
    else
    {
      penalty_upper_bound_constraint = 1;
      penalty_volume_constraint = 1;
      penalty_dtmr_constraint = 1;
    }
  }
  else
  {
    penalty_upper_bound_constraint = 1;
    penalty_volume_constraint = 1;
    penalty_dtmr_constraint = 1;
  }
  total_violation = penalty_upper_bound_constraint * upper_bound_constraint_violation +
                    penalty_volume_constraint * volume_constraint_violation +
                    penalty_dtmr_constraint * dtmr_constraint_violation;

  if(total_violation > brachy_tolerance )
    *constraint_value = total_violation;
  else
    *constraint_value = 0;
  /*printf("upperbound: %lf\n", upper_bound_constraint_violation);
  printf("dtmr: %lf\n", dtmr_constraint_violation);
  printf("volume: %lf\n", volume_constraint_violation);
  for(i = 3; i < number_of_organ_types; i++ )
  {
    for(j = 0; j < number_of_points_in_organ[i]; j++)
    {
      point_id = list_of_dose_calc_points_in_organ[i][j];
      printf("%d: %lf\n",point_id, dose_at_point[point_id]);
    }
  }
  exit(1);*/
  free( dose_at_point );

  if(use_rounding)
  {
    if(revert_to_non_rounding)
    {
      for(i = 0; i < number_of_dwell_positions; i++)
        parameters[i] = backup[i];
    }
    free(backup);
  }
}

void fullEvaluateBrachytherapyTreatmentPlan( double *parameters, double *objective_value, double *constraint_value )
{
  double *dose_at_point, total_violation, obj, volume;
  int     i, organ_id;
  double *backup = NULL;
  double  upper_bound_constraint_violation,
          volume_constraint_violation,
          dtmr_constraint_violation;

  if(use_rounding)
  {
    backup = (double*)Malloc(number_of_dwell_positions*sizeof(double));
    for( i = 0; i < number_of_dwell_positions; i++ )
    {
      backup[i] = parameters[i];
      parameters[i] = (long int)( parameters[i] * 10 + 0.5 ) / 10.0;
    }
  }

  dose_at_point = matrixVectorMultiplication( full_dose_rate_matrix, parameters , full_number_of_dose_calculation_points, number_of_dwell_positions );

  if ( dwell_time_modulation_restriction )
    dtmr_constraint_violation = getDTMRTotalViolation( parameters, 0, 0 );
  else
    dtmr_constraint_violation = 0.0;

  obj = 0.0;
  volume_constraint_violation = 0.0;
  upper_bound_constraint_violation = 0.0;

  for( organ_id = 0; organ_id < number_of_organ_types; organ_id++)
  {
    if( is_treatment_volume[organ_id] )
    {
      obj += getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose(organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0 , 0 );
    }
    else
    {
      volume = getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, lower_bound_for_dose_at_organ[organ_id], 0, 0 );

      if( volume < volume_criteria[organ_id] && fabs(volume - volume_criteria[organ_id] ) >  brachy_tolerance )
        volume_constraint_violation += ( volume_criteria[organ_id] - volume);

      if( use_dmax_constraints )
        upper_bound_constraint_violation += getMaximumDoseInAnOrganTotalViolation( organ_id, full_number_of_points_in_organ[organ_id], full_list_of_dose_calc_points_in_organ[organ_id], dose_at_point, upper_bound_for_dose_at_organ[organ_id]);
    }
  }

  // Convert maximization of V100_PTV to minimization
  *objective_value = -obj;

  // Process constraint violations
  total_violation = upper_bound_constraint_violation +
                    volume_constraint_violation +
                    dtmr_constraint_violation;

  if(total_violation > brachy_tolerance )
    *constraint_value = total_violation;
  else
    *constraint_value = 0;

  free( dose_at_point );

  if(use_rounding)
  {
    if(revert_to_non_rounding)
    {
      for(i = 0; i < number_of_dwell_positions; i++)
        parameters[i] = backup[i];
    }
    free(backup);
  }
}

double getMaximumDoseInAnOrgan( int organ_id, int number_of_dose_calc_points_in_this_organ, int *list_of_dose_calc_points_in_this_organ, double *dose_at_point )
{
  int i, point_id;
  double max_dose = -1;
  for(i = 0; i < number_of_points_in_organ[organ_id]; i++ )
  {
    point_id = list_of_dose_calc_points_in_this_organ[i];
    if(dose_at_point[point_id] > max_dose)
      max_dose = dose_at_point[point_id];
  }

  return max_dose;
}

double getMaximumDoseInAnOrganTotalViolation( int organ_id, int number_of_dose_calc_points_in_this_organ, int *list_of_dose_calc_points_in_this_organ, double *dose_at_point, double maximum_dose )
{
  int i, point_id;
  double total = 0.0;
  for(i = 0; i < number_of_dose_calc_points_in_this_organ; i++)
  {
    point_id = list_of_dose_calc_points_in_this_organ[i];
    if (( dose_at_point[point_id] <= maximum_dose ) || (fabs(dose_at_point[point_id] - maximum_dose) <= brachy_tolerance))
      continue;
    total += (dose_at_point[point_id] - maximum_dose);
  }
  return total;
}

double getRelativeVolumeOfAnOrganReceivingMoreThanOrEqualToACertainDose( int organ_id, int number_of_dose_calc_points_in_this_organ, int *list_of_dose_calc_points_in_this_organ, double *dose_at_point, double minimum_dose, short use_smoothening, double scale_factor )
{
  int i, point_id;
  double relative_volume = 0.0;

  for(i = 0; i < number_of_dose_calc_points_in_this_organ; i++)
  {
    point_id = list_of_dose_calc_points_in_this_organ[i];
    if( (dose_at_point[point_id] >= minimum_dose) || (fabs(dose_at_point[point_id]-minimum_dose) <= brachy_tolerance) )
      relative_volume += 1;
    else
    {
      if(use_smoothening == SMOOTHENING_WITHOUT_SCALE_FACTOR )
        relative_volume += (( 1.0 / ( exp( fabs( dose_at_point[point_id] - minimum_dose ) ) ) )) /
        ( number_of_dose_calc_points_in_this_organ );
      else if(use_smoothening == SMOOTHENING_WITH_SCALE_FACTOR )
      {
        if( scale_factor == 0 )
          scale_factor = 1;
        relative_volume += ( 1.0 / ( scale_factor * exp( fabs( dose_at_point[point_id] - minimum_dose))));
      }
    }
  }

  relative_volume = relative_volume / ((double)number_of_dose_calc_points_in_this_organ);
  return relative_volume;
}

double getRelativeVolumeOfAnOrganReceivingLessThanOrEqualToACertainDose( int organ_id, int number_of_dose_calc_points_in_this_organ, int *list_of_dose_calc_points_in_this_organ, double *dose_at_point, double maximum_dose, short use_smoothening, double scale_factor )
{
  int i, point_id;
  double relative_volume = 0.0;

  for(i = 0; i < number_of_dose_calc_points_in_this_organ; i++)
  {
    point_id = list_of_dose_calc_points_in_this_organ[i];
    if( (dose_at_point[point_id] <= maximum_dose) || (fabs(dose_at_point[point_id]-maximum_dose) <= brachy_tolerance) )
      relative_volume += 1;
    else
    {
      if(use_smoothening == SMOOTHENING_WITHOUT_SCALE_FACTOR )
        relative_volume += (( 1.0 / ( exp( fabs( dose_at_point[point_id] - maximum_dose ) ) ) )) /
        ( number_of_dose_calc_points_in_this_organ );
      else if(use_smoothening == SMOOTHENING_WITH_SCALE_FACTOR )
      {
        if(scale_factor == 0)
          scale_factor = 1;
        relative_volume += ( 1.0 / ( scale_factor * exp( fabs( dose_at_point[point_id] - maximum_dose))));
      }
    }
  }

  relative_volume = relative_volume / ((double)number_of_dose_calc_points_in_this_organ);
  return relative_volume;
}

void expandRandomSetOfDoseCalculationPoints( int *number_of_dose_calc_points )
{
  int i, j, count, organ_id, dwell_id, calc_point_id;

  if( (*number_of_dose_calc_points) > full_number_of_dose_calculation_points )
  {
    printf("The number of requested point is more than the available amount!\n");
    exit(1);
  }

  if( full_list_of_random_order_dose_calc_points_in_organ == NULL )
  {
    full_list_of_random_order_dose_calc_points_in_organ = (int**)Malloc(number_of_organ_types*sizeof(int*));
    for(organ_id  = 0; organ_id < number_of_organ_types; organ_id++)
      full_list_of_random_order_dose_calc_points_in_organ[organ_id] = randomPermutation(full_number_of_points_in_organ[organ_id]);
  }

  ezilaitiniSetOfDoseCalculationPoints();
  number_of_points_in_organ = (int*)Malloc( number_of_organ_types*sizeof(int));
  number_of_dose_calculation_points = 0;
  for(organ_id = 0; organ_id < number_of_organ_types; organ_id++)
  {
    number_of_points_in_organ[organ_id] = (int)round((percentage_of_dose_point_in_organ[organ_id]*(*number_of_dose_calc_points) ));
    number_of_dose_calculation_points += number_of_points_in_organ[organ_id];
    //printf("%d : %d\n", organ_id, number_of_points_in_organ[organ_id]);
  }
  *number_of_dose_calc_points = number_of_dose_calculation_points;
  //printf("Total: %d\n", number_of_dose_calculation_points );
  dose_rate_matrix = (double**)Malloc(number_of_dose_calculation_points*sizeof(double*));
  for(i = 0; i < number_of_dose_calculation_points; i++)
    dose_rate_matrix[i] = (double*)Malloc(number_of_dwell_positions*sizeof(double));

  list_of_dose_calc_points_in_organ = (int**)Malloc(number_of_organ_types*sizeof(int*));
  for(i = 0; i < number_of_organ_types; i++ )
    list_of_dose_calc_points_in_organ[i] = (int*)Malloc(number_of_points_in_organ[i]*sizeof(int));

  count = 0;
  for(organ_id = 0; organ_id < number_of_organ_types; organ_id++)
  {
    for(j = 0; j < number_of_points_in_organ[organ_id]; j++)
    {
      calc_point_id = full_list_of_dose_calc_points_in_organ[organ_id][full_list_of_random_order_dose_calc_points_in_organ[organ_id][j]];
      for(dwell_id = 0; dwell_id < number_of_dwell_positions; dwell_id++)
        dose_rate_matrix[count][dwell_id] = full_dose_rate_matrix[calc_point_id][dwell_id];
      list_of_dose_calc_points_in_organ[organ_id][j] = count;
      count++;
    }
  }
}

void loadRandomSetOfDoseCalculationPoints( int *number_of_dose_calc_points )
{
  int i, j, count, calc_point_index, organ_id, dwell_id;
  int *aList;

  if( (*number_of_dose_calc_points) > full_number_of_dose_calculation_points )
  {
    printf("The number of requested point is more than the available amount!\n");
    exit(1);
  }

  // Free memory for existing set of dose calculation points.
  ezilaitiniSetOfDoseCalculationPoints();

  number_of_points_in_organ = (int*)Malloc( number_of_organ_types*sizeof(int));
  number_of_dose_calculation_points = 0;
  for(organ_id = 0; organ_id < number_of_organ_types; organ_id++)
  {
    number_of_points_in_organ[organ_id] = (int)round((percentage_of_dose_point_in_organ[organ_id]*(*number_of_dose_calc_points) ));
    number_of_dose_calculation_points += number_of_points_in_organ[organ_id];
  }
  *number_of_dose_calc_points = number_of_dose_calculation_points;
  dose_rate_matrix = (double**)Malloc(number_of_dose_calculation_points*sizeof(double*));
  for(i = 0; i < number_of_dose_calculation_points; i++)
    dose_rate_matrix[i] = (double*)Malloc(number_of_dwell_positions*sizeof(double));

  list_of_dose_calc_points_in_organ = (int**)Malloc(number_of_organ_types*sizeof(int*));
  for(i = 0; i < number_of_organ_types; i++ )
    list_of_dose_calc_points_in_organ[i] = (int*)Malloc(number_of_points_in_organ[i]*sizeof(int));

  count = 0;
  for(organ_id = 0; organ_id < number_of_organ_types; organ_id++)
  {
    aList = getRandomSetOfDoseCalculationPointsFromAnOrgan( organ_id, number_of_points_in_organ[organ_id]);
    for(j = 0; j < number_of_points_in_organ[organ_id]; j++)
    {
      calc_point_index = aList[j];
      for(dwell_id = 0; dwell_id < number_of_dwell_positions; dwell_id++)
        dose_rate_matrix[count][dwell_id] = full_dose_rate_matrix[calc_point_index][dwell_id];
      list_of_dose_calc_points_in_organ[organ_id][j] = count;
      count++;
    }
    free(aList);
  }

  FILE *fp;
  fp = fopen("./dose_rate_matrix.txt", "w");
  for( i = 0; i < number_of_dose_calculation_points; i++ )
  {
    for( j = 0; j < number_of_dwell_positions; j++ )
      fprintf( fp, "%lf ", dose_rate_matrix[i][j] );
    fprintf( fp, "\n" );
  }
  fclose( fp );
}

int* getRandomSetOfDoseCalculationPointsFromAnOrgan( int organ_id, int number_of_points )
{
  int *result, *order, i;

  if( number_of_points > full_number_of_points_in_organ[organ_id] )
  {
    printf("The number of requested point is more than the available amount!\n");
    exit(1);
  }
  result = (int*) Malloc(number_of_points*sizeof(int));
  order = randomPermutation(full_number_of_points_in_organ[organ_id]);
  for( i  = 0; i < number_of_points; i++ )
    result[i] = full_list_of_dose_calc_points_in_organ[organ_id][order[i]];

  free( order );
  return result;
}

double getDTMRTotalViolation( double *parameters, short use_dynamic_evaluation, double scale_factor )
{
  int i, j, index, index_next;
  double time1, time2, max_time, min_time;
  double total, g;
  total = 0;
  for( i = 0; i < number_of_catheters; i++ )
  {
    for( j = 0; j < number_of_dwell_positions_in_catheter[i]-1; j++ )
    {
      index = list_of_dwell_positions_in_catheter[i][j];
      index_next = list_of_dwell_positions_in_catheter[i][j+1];
      time1 = parameters[index];
      time2 = parameters[index_next];
      max_time = max( time1, time2 );
      min_time = min( time1, time2 );

      if( !use_dynamic_evaluation )
      {
        if( (max_time <= 2.0 *min_time + brachy_tolerance ) || (max_time - min_time <= 0.1 + brachy_tolerance ))
          continue;
        total += ( max_time - 2.0 *min_time);
      }
      else
      {
        if( scale_factor > 0 )
          g = scale_factor;
        else
          g = 1;
        if( (max_time <= (2.0 + 100.0/(g*g))*min_time + brachy_tolerance ) || (max_time - min_time <= 0.1 + brachy_tolerance ))
          continue;
        total += ( max_time - (2.0 + 100.0/(g*g))*min_time );
      }
    }
  }
  return total;
}

void readReferenceFileForDVHBasedInversePlanning( char *filename , char *distance_filename, int *number_of_parameters )
{
  FILE *file;
  int i, j, type_of_organ, index, catheter_index, pos_index, temp_int, number_of_points;
  int maximum_number_of_dwell_positions_in_a_catheter;
  double temp;

  file = NULL;
  file = fopen( filename, "r" );
  if ( file == NULL )
  {
    printf( "Cannot read input file!\n " );
    exit(1);
  }

  fscanf( file, "%d", &number_of_catheters);                // read the number of catheters
  fscanf( file, "%d", &number_of_dwell_positions );         // read the number of dwell positions;

  *number_of_parameters = number_of_dwell_positions;
  // read each dwell position and its corresponding catheter
  number_of_dwell_positions_in_catheter = ( int* ) Malloc( number_of_catheters * sizeof( int ) );
  for(i = 0; i < number_of_catheters; i++ )
    number_of_dwell_positions_in_catheter[i] = 0;
  catheter_of_dwell_position = ( int* ) Malloc( number_of_dwell_positions * sizeof( int ) );
  list_of_dwell_positions_in_catheter = ( int** ) Malloc( number_of_catheters * sizeof( int* ) );
  maximum_number_of_dwell_positions_in_a_catheter = 100;
  for( i = 0; i < number_of_catheters; i++ )
  {
    list_of_dwell_positions_in_catheter[i] = ( int* ) Malloc( maximum_number_of_dwell_positions_in_a_catheter * sizeof( int ) );
    for( j = 0; j < maximum_number_of_dwell_positions_in_a_catheter; j++ )
      list_of_dwell_positions_in_catheter[i][j] = -1;
  }

  for( i = 0; i < number_of_dwell_positions; i++ )
  {
    fscanf( file, "%d", &pos_index );
    pos_index = pos_index - 1;
    fscanf( file, "%d", &catheter_index );
    catheter_index = catheter_index - 1;

    catheter_of_dwell_position[pos_index] = catheter_index;
    list_of_dwell_positions_in_catheter[catheter_index][number_of_dwell_positions_in_catheter[catheter_index]] = pos_index;
    number_of_dwell_positions_in_catheter[catheter_index] += 1;
  }

  fscanf( file, "%d", &full_number_of_dose_calculation_points ); // read the number of dose calculation points;

  // read the dose-rate matrix: nrow = number_of_dose_calculation_points; ncol = number_of_dwell_positions
  full_dose_rate_matrix = ( double** ) Malloc( full_number_of_dose_calculation_points * sizeof ( double* ) );
  for( i = 0; i < full_number_of_dose_calculation_points; i++ )
    full_dose_rate_matrix[i] = ( double* ) Malloc( number_of_dwell_positions * sizeof( double ) );

  for( i = 0; i < full_number_of_dose_calculation_points; i++ )
    for( j = 0; j < number_of_dwell_positions; j++ )
    {
      fscanf( file, "%lf", &temp);
      full_dose_rate_matrix[i][j] =  temp/100.0;
    }

  // read the number of organ types considered during inverse planning
  fscanf( file, "%d", &number_of_organ_types );
  full_number_of_points_in_organ = ( int* ) Malloc( number_of_organ_types * sizeof( int ) );
  for( i = 0; i < number_of_organ_types; i++ )
    full_number_of_points_in_organ[i] = 0;

  // read the type of organ and the list of dose calc points in that organ
  percentage_of_dose_point_in_organ = (double*)Malloc(number_of_organ_types*sizeof(double));
  full_list_of_dose_calc_points_in_organ = ( int** ) Malloc( number_of_organ_types * sizeof( int* ));
  for( i = 0; i < number_of_organ_types; i++ )
  {
    fscanf( file, "%d", &type_of_organ );
    fscanf( file, "%d", &number_of_points );
    type_of_organ = type_of_organ - 1; // C uses index-0 while data files use index-1
    full_number_of_points_in_organ[type_of_organ] = number_of_points;
    full_list_of_dose_calc_points_in_organ[type_of_organ] = (int*) Malloc(number_of_points*sizeof(int));
    for( j = 0; j < number_of_points; j++ )
    {
      fscanf( file, "%d", &index );
      full_list_of_dose_calc_points_in_organ[type_of_organ][j] = index - 1;
    }
    percentage_of_dose_point_in_organ[type_of_organ] = ((double)number_of_points) / ((double)full_number_of_dose_calculation_points);
    //printf("%d : %lf - %d\n", type_of_organ, percentage_of_dose_point_in_organ[type_of_organ], number_of_points);
  }

  organ_volume                  = ( double* ) Malloc( number_of_organ_types * sizeof( double ) );
  for( i = 0; i < number_of_organ_types; i++ )
  {
    fscanf( file, "%d", &temp_int);
    fscanf( file, "%lf", &temp);
    temp_int = temp_int - 1;
    organ_volume[temp_int] = temp; // volume of the organ in cc
  }
  // read dvh index
  fscanf( file, "%d", &number_of_dvh_indices );
  fscanf( file, "%lf", &temp);
  prescribed_dose = temp;
  dvh_index_dose_level          = ( double* ) Malloc( number_of_dvh_indices * sizeof( double ) );
  dvh_index_upper_bound         = ( double* ) Malloc( number_of_dvh_indices * sizeof( double ) );
  organ_id_of_dvh_index         = ( int* ) Malloc (number_of_dvh_indices * sizeof( int ) );
  is_treatment_index            = ( int*) Malloc( number_of_dvh_indices * sizeof( int ) );
  is_dose_index                 = ( int*) Malloc( number_of_dvh_indices * sizeof( int ) );

  for( i = 0; i < number_of_dvh_indices; i++ )
  {
    fscanf( file, "%d", &type_of_organ );
    type_of_organ = type_of_organ - 1;
    organ_id_of_dvh_index[i] = type_of_organ;
    fscanf( file, "%lf", &temp );
    dvh_index_dose_level[i] = (temp/100.0)*prescribed_dose;
    fscanf( file, "%lf", &temp );
    dvh_index_upper_bound[i] = temp;
    fscanf( file, "%d", &temp_int );
    is_dose_index[i] = temp_int;
    fscanf( file, "%d", &temp_int );
    is_treatment_index[i] = temp_int;
    if( is_treatment_index[i] )
      treatment_volume_id = type_of_organ;
    if( is_dose_index[i] )
    {
      dvh_index_upper_bound[i] = dvh_index_upper_bound[i] / organ_volume[type_of_organ];
      //printf("index: %d  dvh_upper_bound: %lf\n",i, dvh_index_upper_bound[i]);
    }
  }
  fclose( file );

  distance_matrix = NULL;
  if( distance_filename == NULL )
    return;

  file = fopen( distance_filename, "r");
  if( file == NULL )
  {
    printf("distance matrix not exit!\n");
    return;
  }
  distance_matrix  = (double**) Malloc(number_of_dwell_positions*sizeof( double* ));
  for( i = 0; i < number_of_dwell_positions; i++ )
  {
    distance_matrix[i] = (double*) Malloc(number_of_dwell_positions*sizeof( double ));
    for( j = 0; j < number_of_dwell_positions; j++ )
    {
      fscanf( file, "%lf", &temp );
      distance_matrix[i][j] = temp;
    }
  }
  fclose( file );
}
void readReferenceFileForBrachytherapyLinearVolumeBasedInversePlanning( char *filename , int *number_of_parameters )
{
  FILE *file;
  int i, j, type_of_organ, index, catheter_index, pos_index, temp_int, number_of_points;
  int maximum_number_of_dwell_positions_in_a_catheter;
  double temp;

  file = NULL;
  file = fopen( filename, "r" );
  if ( file == NULL )
  {
    printf( "Cannot read input file!\n " );
    exit(1);
  }

  fscanf( file, "%d", &number_of_catheters);                // read the number of catheters
  fscanf( file, "%d", &number_of_dwell_positions );         // read the number of dwell positions;

  *number_of_parameters = number_of_dwell_positions;
  // read each dwell position and its corresponding catheter
  number_of_dwell_positions_in_catheter = ( int* ) Malloc( number_of_catheters * sizeof( int ) );
  for(i = 0; i < number_of_catheters; i++ )
    number_of_dwell_positions_in_catheter[i] = 0;
  catheter_of_dwell_position = ( int* ) Malloc( number_of_dwell_positions * sizeof( int ) );
  list_of_dwell_positions_in_catheter = ( int** ) Malloc( number_of_catheters * sizeof( int* ) );
  maximum_number_of_dwell_positions_in_a_catheter = 100;
  for( i = 0; i < number_of_catheters; i++ )
  {
    list_of_dwell_positions_in_catheter[i] = ( int* ) Malloc( maximum_number_of_dwell_positions_in_a_catheter * sizeof( int ) );
    for( j = 0; j < maximum_number_of_dwell_positions_in_a_catheter; j++ )
      list_of_dwell_positions_in_catheter[i][j] = -1;
  }

  for( i = 0; i < number_of_dwell_positions; i++ )
  {
    fscanf( file, "%d", &pos_index );
    pos_index = pos_index - 1;
    fscanf( file, "%d", &catheter_index );
    catheter_index = catheter_index - 1;

    catheter_of_dwell_position[pos_index] = catheter_index;
    list_of_dwell_positions_in_catheter[catheter_index][number_of_dwell_positions_in_catheter[catheter_index]] = pos_index;
    number_of_dwell_positions_in_catheter[catheter_index] += 1;
  }

  fscanf( file, "%d", &full_number_of_dose_calculation_points ); // read the number of dose calculation points;

  // read the dose-rate matrix: nrow = number_of_dose_calculation_points; ncol = number_of_dwell_positions
  full_dose_rate_matrix = ( double** ) Malloc( full_number_of_dose_calculation_points * sizeof ( double* ) );
  for( i = 0; i < full_number_of_dose_calculation_points; i++ )
    full_dose_rate_matrix[i] = ( double* ) Malloc( number_of_dwell_positions * sizeof( double ) );

  for( i = 0; i < full_number_of_dose_calculation_points; i++ )
    for( j = 0; j < number_of_dwell_positions; j++ )
    {
      fscanf( file, "%lf", &temp);
      full_dose_rate_matrix[i][j] =  temp/100.0;
    }

  // read the number of organ types considered during inverse planning
  fscanf( file, "%d", &number_of_organ_types );
  full_number_of_points_in_organ = ( int* ) Malloc( number_of_organ_types * sizeof( int ) );
  for( i = 0; i < number_of_organ_types; i++ )
    full_number_of_points_in_organ[i] = 0;

  // read the type of organ and the list of dose calc points in that organ
  percentage_of_dose_point_in_organ = (double*)Malloc(number_of_organ_types*sizeof(double));
  full_list_of_dose_calc_points_in_organ = ( int** ) Malloc( number_of_organ_types * sizeof( int* ));
  for( i = 0; i < number_of_organ_types; i++ )
  {
    fscanf( file, "%d", &type_of_organ );
    fscanf( file, "%d", &number_of_points );
    type_of_organ = type_of_organ - 1; // C uses index-0 while data files use index-1
    full_number_of_points_in_organ[type_of_organ] = number_of_points;
    full_list_of_dose_calc_points_in_organ[type_of_organ] = (int*) Malloc(number_of_points*sizeof(int));
    for( j = 0; j < number_of_points; j++ )
    {
      fscanf( file, "%d", &index );
      full_list_of_dose_calc_points_in_organ[type_of_organ][j] = index - 1;
    }
    percentage_of_dose_point_in_organ[type_of_organ] = ((double)number_of_points) / ((double)full_number_of_dose_calculation_points);
    //printf("%d : %lf - %d\n", type_of_organ, percentage_of_dose_point_in_organ[type_of_organ], number_of_points);
  }

  lower_bound_for_dose_at_organ =       ( double* ) Malloc( number_of_organ_types * sizeof( double ) );
  upper_bound_for_dose_at_organ =       ( double* ) Malloc( number_of_organ_types * sizeof( double ) );
  volume_criteria               =       ( double* ) Malloc( number_of_organ_types * sizeof( double ) );
  is_treatment_volume           =       ( int* ) Malloc( number_of_organ_types * sizeof( int ) );

  for( i = 0; i < number_of_organ_types; i++ )
  {
    lower_bound_for_dose_at_organ[i] = 0;
    upper_bound_for_dose_at_organ[i] = 1e+308;
    volume_criteria[i] = 0;
  }

  for( i = 0; i < number_of_organ_types; i++ )
  {
    fscanf( file, "%d", &type_of_organ );
    type_of_organ = type_of_organ - 1;

    fscanf( file, "%lf", &temp );
    if( temp >= 0 )
      lower_bound_for_dose_at_organ[type_of_organ] = temp;
    else
      lower_bound_for_dose_at_organ[type_of_organ] = 0;

    fscanf( file, "%lf", &temp );
    if( temp >= 0 )
      upper_bound_for_dose_at_organ[type_of_organ] = temp;
    else
      upper_bound_for_dose_at_organ[type_of_organ] = 1e+8;

    fscanf( file, "%lf", &temp );
    if( temp >= 0 )
      volume_criteria[type_of_organ] = temp;
    else
      volume_criteria[type_of_organ] = 0;

    fscanf( file, "%d", &temp_int);
    is_treatment_volume[type_of_organ] = temp_int;
    if( is_treatment_volume[type_of_organ] )
      treatment_volume_id = type_of_organ;
  }
  fclose( file );
}

void ezilaitiniSetOfDoseCalculationPoints()
{
  int i;
  if( dose_rate_matrix == NULL )
    return;

  for( i = 0; i < number_of_dose_calculation_points; i++ )
    free( dose_rate_matrix[i] );
  free( dose_rate_matrix );

  for( i = 0; i < number_of_organ_types; i++ )
    free( list_of_dose_calc_points_in_organ[i] );
  free( list_of_dose_calc_points_in_organ );

  free( number_of_points_in_organ );
}
void ezilaitiniBrachyTherapyDVH()
{
  int i;

  for( i = 0; i < number_of_dose_calculation_points; i++ )
    free( dose_rate_matrix[i] );
  free( dose_rate_matrix );

  for( i = 0; i < number_of_organ_types; i++ )
    free( list_of_dose_calc_points_in_organ[i] );
  free( list_of_dose_calc_points_in_organ );

  free( number_of_points_in_organ );

  free( is_treatment_index );
  free( is_dose_index );
  free( dvh_index_dose_level );
  free( dvh_index_upper_bound );
  free( organ_id_of_dvh_index );
  free( organ_volume );

  free( catheter_of_dwell_position );
  free( number_of_dwell_positions_in_catheter );
  for( i = 0; i < number_of_catheters; i++ )
    free( list_of_dwell_positions_in_catheter[i] );
  free( list_of_dwell_positions_in_catheter );

  for( i = 0; i < full_number_of_dose_calculation_points; i++ )
    free( full_dose_rate_matrix[i] );
  free( full_dose_rate_matrix );

  for( i = 0; i < number_of_organ_types; i++ )
    free( full_list_of_dose_calc_points_in_organ[i] );
  free( full_list_of_dose_calc_points_in_organ );

  free( full_number_of_points_in_organ );
  free( percentage_of_dose_point_in_organ );
  if( full_list_of_random_order_dose_calc_points_in_organ != NULL )
  {
    for( i = 0; i < number_of_organ_types; i++ )
      free( full_list_of_random_order_dose_calc_points_in_organ[i] );
    free( full_list_of_random_order_dose_calc_points_in_organ );
  }

  if( distance_matrix != NULL )
  {
    for( i = 0; i < number_of_dwell_positions; i++ )
      free( distance_matrix[i] );
    free( distance_matrix );
  }
}


void ezilaitiniBrachyTherapy()
{
  int i;

  for( i = 0; i < number_of_dose_calculation_points; i++ )
    free( dose_rate_matrix[i] );
  free( dose_rate_matrix );

  for( i = 0; i < number_of_organ_types; i++ )
    free( list_of_dose_calc_points_in_organ[i] );
  free( list_of_dose_calc_points_in_organ );

  free( number_of_points_in_organ );

  free( lower_bound_for_dose_at_organ );
  free( upper_bound_for_dose_at_organ );
  free( volume_criteria );

  free( catheter_of_dwell_position );
  free( number_of_dwell_positions_in_catheter );
  for( i = 0; i < number_of_catheters; i++ )
    free( list_of_dwell_positions_in_catheter[i] );
  free( list_of_dwell_positions_in_catheter );

  free( is_treatment_volume );

  for( i = 0; i < full_number_of_dose_calculation_points; i++ )
    free( full_dose_rate_matrix[i] );
  free( full_dose_rate_matrix );

  for( i = 0; i < number_of_organ_types; i++ )
    free( full_list_of_dose_calc_points_in_organ[i] );
  free( full_list_of_dose_calc_points_in_organ );

  free( full_number_of_points_in_organ );
  free( percentage_of_dose_point_in_organ );
  if( full_list_of_random_order_dose_calc_points_in_organ != NULL )
  {
    for( i = 0; i < number_of_organ_types; i++ )
      free( full_list_of_random_order_dose_calc_points_in_organ[i] );
    free( full_list_of_random_order_dose_calc_points_in_organ );
  }
}
