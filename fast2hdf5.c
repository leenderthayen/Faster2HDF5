//-------------------------------------------------------//
//     ADC_MOSAHR Data converstion from .fast to .h5
//-------------------------------------------------------//

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <sys/stat.h>

#include "fasterac/fasterac.h"
#include "fasterac/fast_data.h"
#include "fasterac/utils.h"
#include "fasterac/adc_caras.h"

#include "hdf5.h"
#include "hdf5_hl.h"

#define NFA_EVT 5 /*number of variables for ADC type event*/
#define NFA_CNT 4 /*number of variables for A_COUNT type event*/
#define BUFSIZE 10000
#define NROFDETECTORS 3

//-------------------------------------------------------//

void display_usage (char* prog_name) {
    printf ("\n");
    printf ("***  %s  *** \n", prog_name);
    printf ("\n");
    printf ("Convert faster data files to HDF5 format.\n");
    printf ("\n");
    printf ("usage : \n");
    printf ("       %s  input_file.fast  output_file.h5 [label_d1 ... label_d4]\n", prog_name);
    printf ("\n");
    printf ("The input_file.fast file is converted to a output_file.h5 file.\n");
    printf ("label_d1 and label_d4 denote the physical channels on the MOSAHR daughterboard for all detectors.\n");
    printf ("Data is split up for all detectors, with separate tables for ADC and Counter data.\n");
    printf ("\n");
    printf ("\n");
}

//-------------------------------------------------------//

int main (int argc, char** argv) {

  //------------------------------//
  //   Variable definition 
  //------------------------------//
  faster_file_reader_p reader;                                //  data file reader
  faster_data_p        data;                                  //  data pointer
  char                 prog_name  [256];
  char                 input_file [256];
  char                 output_file[256];
  sampling             s;
  adc_data             a;
  adc_counter          a_count;
  int                  indices[2*NROFDETECTORS];
  //  faster data  (fasterac.h)
  unsigned char      alias;
  unsigned short     label;
  unsigned short     lsize;
  unsigned long long clock;           // time stamp ns
  long double        hr_clock;        // time stamp + tdc
  unsigned short labels[NROFDETECTORS];

  hid_t              file;                          /* handles */
  hid_t		     groups[NROFDETECTORS];
  hid_t              at[NROFDETECTORS];
  hid_t              sds;
  hsize_t    chunk_size = 40;
  int        *fill_data = NULL;
  int        compress  = 0;
  herr_t     status;
  int        i;


  strcpy (prog_name,  argv[0]);                                //  prog_name

  if (argc < 4) {                                              //  command args & usage
    printf("Not enough arguments\n");
    display_usage (prog_name);
    return EXIT_FAILURE;
  }

  else {
    for (int i = 3 ; i < argc; i++) {
      if (i-3 < NROFDETECTORS) {
        labels[i-3] = atoi (argv[i]);
      }
    }
  }

  strcpy (input_file, argv[1]);                                //  input file
  strcpy (output_file, argv[2]);

  reader = faster_file_reader_open (input_file);               //  create a reader
  if (reader == NULL) {
    printf ("ERROR opening file %s\n", input_file);
    display_usage (prog_name);
    return EXIT_FAILURE;
  }

  if (labels[0] == 0){
    printf("Please enter the channel for Detector 1: ");
    scanf("%hu", &labels[0]);
  }

 //TODO: struct for OSC data


 //----------------------------------------------//
 //     Table definitions for ADC type event
 //----------------------------------------------//

 /* Definition of variables */
 typedef struct a_evt
 {
  signed long long   clock;
  int    delta_t;
  int    adc;
  unsigned int pileup;
  unsigned short saturation;
 } a_evt;

 /* Calculate the size and the offsets of our struct members in memory */
 a_evt dst_buf_evt[1];

 size_t dst_size_evt =  sizeof( a_evt );
 size_t dst_offset_evt[NFA_EVT] = {HOFFSET(a_evt, clock),
				HOFFSET(a_evt, delta_t),
				HOFFSET(a_evt, adc),
				HOFFSET(a_evt, pileup),
				HOFFSET(a_evt, saturation)};

 size_t dst_sizes_evt[NFA_EVT] = {sizeof( dst_buf_evt[0].clock),
                               sizeof( dst_buf_evt[0].delta_t),
                               sizeof( dst_buf_evt[0].adc),
			       sizeof( dst_buf_evt[0].pileup),
			       sizeof( dst_buf_evt[0].saturation)};

  /* Define field information */
  const char *field_names_evt[NFA_EVT]  =
  { "Clock", "Delta_t", "ADC", "Pileup", "Saturation" };

  hid_t      field_type_evt[NFA_EVT];

  /* Initialize field_type */
  field_type_evt[0] = H5T_NATIVE_LLONG;
  field_type_evt[1] = H5T_NATIVE_UINT;
  field_type_evt[2] = H5T_NATIVE_UINT;
  field_type_evt[3] = H5T_NATIVE_USHORT;
  field_type_evt[4] = H5T_NATIVE_USHORT;


 //----------------------------------------------//
 //     Table definitions for A_COUNT type event
 //----------------------------------------------//

 /* Definition of variables */
 typedef struct
 {
  signed long long clock;
  int calc;
  int sent;
  int trig;
 } a_cnt;

 /* Calculate the size and the offsets of our struct members in memory */
 a_cnt dst_buf_cnt[1];

 size_t dst_size_cnt =  sizeof( a_cnt );
 size_t dst_offset_cnt[NFA_CNT] = {HOFFSET(a_cnt, clock),
				HOFFSET(a_cnt, calc),
				HOFFSET(a_cnt, sent),
				HOFFSET(a_cnt, trig)};

 size_t dst_sizes_cnt[NFA_CNT] = {sizeof( dst_buf_cnt[0].clock),
                               sizeof( dst_buf_cnt[0].calc),
                               sizeof( dst_buf_cnt[0].sent),
			       sizeof( dst_buf_cnt[0].trig)};

  /* Define field information */
  const char *field_names_cnt[NFA_CNT]  =
  {"Clock", "Calc", "Sent", "Trig"};

  hid_t      field_type_cnt[NFA_CNT];

  /* Initialize field_type */
  field_type_cnt[0] = H5T_NATIVE_LLONG;
  field_type_cnt[1] = H5T_NATIVE_UINT;
  field_type_cnt[2] = H5T_NATIVE_UINT;
  field_type_cnt[3] = H5T_NATIVE_UINT;

  /* Save old error handler */
  herr_t (*old_func)(void*);
  void *old_client_data;
  H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
  /* Turn off error handling */
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  file = H5Fopen(output_file, H5F_ACC_RDWR, H5P_DEFAULT);

  /* Restore previous error handler */
  H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);

  if (file < 0){
    printf("File does not exist.\n");

    /*
     * Create a new file. If file exists its contents will be overwritten.
     */
    file = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // attribute dataspace
    hsize_t x = 1;
    sds = H5Screate_simple(1, &x, NULL);

    char groupname[20];
    char tablenameADC[30];
    char tablenameCounters[30];

    // Group Creation
    for (int i = 0; i < NROFDETECTORS; ++i) {
      if (labels[i] != 0) {
        sprintf(groupname, "Detector%d", (i+1));
        // Create group
        groups[i] = H5Gcreate2(file, groupname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        at[i] = H5Acreate(groups[i], "Channel", H5T_NATIVE_USHORT, sds, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite(at[i], H5T_NATIVE_USHORT, &labels[i]);
        // Create ADC table
        sprintf(tablenameADC, "Detector %d ADC", (i+1));
        status = H5TBmake_table( tablenameADC, groups[i], "ADC", NFA_EVT, 0,
                         dst_size_evt, field_names_evt, dst_offset_evt, field_type_evt,
                         chunk_size, fill_data, compress, NULL);
        // Create Counters table
        sprintf(tablenameCounters, "Detector %d Counters", (i+1));
        status=H5TBmake_table( tablenameCounters, groups[i], "Counters", NFA_CNT, 0,
                         dst_size_cnt, field_names_cnt, dst_offset_cnt, field_type_cnt,
                         chunk_size, fill_data, compress, NULL);
      }
    }
  }
  else {
    char groupname[20];
    for (int i = 0; i < NROFDETECTORS; ++i) {
      if (labels[i] != 0) {
        sprintf(groupname, "Detector%d", (i+1));
        groups[i] = H5Gopen(file, groupname, H5P_DEFAULT);
        at[i] = H5Aopen(groups[i], "Channel", H5P_DEFAULT);
      }
    }
  }


  a_evt buf_evt[NROFDETECTORS][BUFSIZE];

  a_cnt buf_cnt[NROFDETECTORS][BUFSIZE];

  int j = 0;

  // Main Loop
  while ((data = faster_file_reader_next (reader)) != NULL) {
    alias = faster_data_type_alias (data);
    label = faster_data_label      (data);
    clock = faster_data_clock_ns   (data);
    lsize = faster_data_load_size  (data);

    for (int i = 0; i < NROFDETECTORS; ++i) {
      if (labels[i] != 0) {
        if (indices[2*i] > 0 && indices[2*i] == BUFSIZE) {
           printf("Cycle %d\n", j);
           status = H5TBappend_records(groups[i], "ADC", BUFSIZE, dst_size_evt, dst_offset_evt, dst_sizes_evt,
                                       &buf_evt[i][0]);
           indices[2*i] = 0;
        }
        if (indices[2*i+1] > 0 && indices[2*i+1] == BUFSIZE){
          printf("Cycle %d\n", j);
          status = H5TBappend_records(groups[i], "Counters", BUFSIZE, dst_size_cnt, dst_offset_cnt,
                                      dst_sizes_cnt, &buf_cnt[i][0]);
          indices[2*i+1] = 0;
        }
      }
    }

    switch (alias) {
       case SAMPLING_TYPE_ALIAS:
          faster_data_load(data, &s);
          //do something
          break;
       case ADC_DATA_TYPE_ALIAS:
          faster_data_load(data, &a);
	  a_evt d = {clock, adc_delta_t_ns(a), a.measure, a.pileup, a.saturated};
          for (int i = 0; i < NROFDETECTORS; ++i) {
            if (labels[i] == label) {
              buf_evt[i][indices[2*i]] = d;
              indices[2*i] = indices[2*i]+1;
            }
          }
          j = j+1;
          break;
       case ADC_COUNTER_TYPE_ALIAS:
          faster_data_load(data, &a_count);
          a_cnt d_cnt = {clock, a_count.calc, a_count.sent, a_count.trig};
          for (int i = 0; i < NROFDETECTORS; ++i) {
            if (labels[i] == label%1000) {
              buf_cnt[i][indices[2*i+1]] = d_cnt;
              indices[2*i+1] = indices[2*i+1]+1;
            }
          }
          j = j+1;
          break;
       default:
         break;
    }
  }

  for (int i = 0; i < NROFDETECTORS; ++i) {
    if (indices[2*i] > 0) {
      a_evt buf[indices[2*i]];
      memcpy(&buf, &buf_evt[i][0], indices[2*i]*sizeof(buf_evt[i][0]));
      status = H5TBappend_records(groups[i], "ADC", indices[2*i], dst_size_evt, dst_offset_evt, dst_sizes_evt,
                                  buf);
    }
    if (indices[2*i+1] > 0) {
      a_cnt buf[indices[2*i+1]];
      memcpy(&buf, &buf_cnt[i][0], indices[2*i+1]*sizeof(buf_cnt[i][0]));
      status = H5TBappend_records(groups[i], "Counters", indices[2*i+1], dst_size_cnt, dst_offset_cnt,
                                  dst_sizes_cnt, buf);
    }
  }

  faster_file_reader_close (reader);                           //  close the reader
  for (int i = 0; i < NROFDETECTORS; ++i) {
    if (labels[i] != 0) {
      H5Aclose(at[i]);
      H5Gclose(groups[i]);
    }
  }
  H5Fclose(file);
  return EXIT_SUCCESS;
}
