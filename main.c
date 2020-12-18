
/*********************************************************************************************

    This is public domain software that was developed by or for the U.S. Naval Oceanographic
    Office and/or the U.S. Army Corps of Engineers.

    This is a work of the U.S. Government. In accordance with 17 USC 105, copyright protection
    is not available for any work of the U.S. Government.

    Neither the United States Government, nor any employees of the United States Government,
    nor the author, makes any warranty, express or implied, without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, or assumes any liability or
    responsibility for the accuracy, completeness, or usefulness of any information,
    apparatus, product, or process disclosed, or represents that its use would not infringe
    privately-owned rights. Reference herein to any specific commercial products, process,
    or service by trade name, trademark, manufacturer, or otherwise, does not necessarily
    constitute or imply its endorsement, recommendation, or favoring by the United States
    Government. The views and opinions of authors expressed herein do not necessarily state
    or reflect those of the United States Government, and shall not be used for advertising
    or product endorsement purposes.

*********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <getopt.h>

#include "nvutility.h"

#include "misp.h"
#include "chrtr2.h"

#include "version.h"


#define         FILTER 9



/*

    Programmer : Jan C. Depner
    Date : 12/28/10

    This program is a filter that uses a low resolution grid (at 5 times cell spacing) of a CHRTR2 file.  First
    we create the low resolution grid at 5 times the spacing of the original grid using every point in the original
    grid as input to MISP.  After we have created the low resolution grid we compare it to the original input grid.
    If all the data in a given cell for the low resolution grid is interpolated, we load the interpolated value
    into MISP.  If there is a real (or digitized) point in the original grid area that corresponds to the
    low resolution grid cell and the depth of the low resolution grid cell is less than 500 fathoms (at 4800 feet per
    second) we load the value into MISP.  If there is a real (or digitized) point in the original grid area
    that corresponds to the low resolution grid cell we check to see if the high and low points in the original
    grid are more than 10 fathoms (at 4800 feet per second) from the low resolution value.  If either exceeeds this
    threshold they are loaded into MISP.  Once all of the low resolution grid cells are tested in this way and all
    of the points are input to MISP we run MISP to produce a new CHRTR2 file at the resolution of the
    original grid.  We read the original input file to retrieve the horizontal uncertainty, vertical uncertainty, 
    status, uncertainty, and number of points (all of these fields may not be set).  We write the new CHRTR2 file
    using all of the above information but replacing the Z value with the newly interpolated (smoothed) Z value.
    We then run through the list of points that were loaded into MISP from the original grid and set the value and
    status of the corresponding cells in the new grid to these values.  That prevents MISP from shaving off the tops
    or bottoms of features.

*/


/*  Structure used for saving load points.  */

typedef struct
{
  NV_I32_COORD2      coord;
  NV_F64_COORD3      xyz;
} LOAD_POINT;



void usage ()
{
  fprintf (stderr, "\n\nUsage: chrtr2_filter5 CHRTR2_FILENAME\n\n");
  fprintf (stderr, "\tThis will produce a filtered CHRTR2 file with a\n");
  fprintf (stderr, "\tname of CHRTR2_FILENAME__filter5.ch2\n\n");
  fflush (stderr);
  exit (-1);
}



void get_load_values (NV_F64_COORD3 grid5_xyz, int32_t row, int32_t col, int32_t chrtr2_handle, CHRTR2_HEADER chrtr2_header, LOAD_POINT *load_point,
                      int32_t *count)
{
  int32_t            i, j, start_row, end_row, start_col, end_col;
  NV_I32_COORD2      coord;
  CHRTR2_RECORD      chrtr2_record;
  NV_F64_COORD3      xyz;
  LOAD_POINT         load_min, load_max;
  uint8_t            real;


  /*  Zero the count.  */

  *count = 0;


  /*  Define the start and end rows and columns for the original CHRTR2 file that correspond to the low resolution cell.  */

  start_row = (row - FILTER) * 5;
  end_row = start_row + 5;
  start_col = (col - FILTER) * 5;
  end_col = start_col + 5;

  memset (&load_min, 0, sizeof (LOAD_POINT));
  memset (&load_max, 0, sizeof (LOAD_POINT));

  load_min.xyz.z = 999999999.0;
  load_max.xyz.z = -999999999.0;
  real = NVFalse;


  /*  Check all of the original grid cells against the low resolution cell.  */

  for (i = start_row ; i < end_row ; i++)
    {
      coord.y = i;

      if (coord.y >= 0 && coord.y < chrtr2_header.height)
        { 
          for (j = start_col ; j < end_col ; j++)
            {
              coord.x = j;

              if (coord.x >= 0 && coord.x < chrtr2_header.width)
                {
                  if (chrtr2_read_record (chrtr2_handle, coord, &chrtr2_record))
                    {
                      chrtr2_perror ();
                      exit (-1);
                    }


                  if (chrtr2_record.status & (CHRTR2_REAL | CHRTR2_DIGITIZED_CONTOUR)) real = NVTrue;


                  chrtr2_get_lat_lon (chrtr2_handle, &xyz.y, &xyz.x, coord);
                  xyz.z = chrtr2_record.z;


                  /*  Get the min and max.  */

                  if (chrtr2_record.z < load_min.xyz.z)
                    {
                      load_min.coord = coord;
                      load_min.xyz = xyz;
                    }

                  if (chrtr2_record.z > load_max.xyz.z)
                    {
                      load_max.coord = coord;
                      load_max.xyz = xyz;
                    }
                }
            }
        }
    }


  /*  If we had real or hand drawn data we have to do the depth checks.  */

  if (real)
    {
      /*  If the average depth is less than 500 fathoms (at 4800 feet per second) or the difference between the average depth
          and the minimum depth is greater than 9.9 fathoms (4800fps) we save the minimum location and depth.  */

      if (grid5_xyz.z * 0.533333 < 500.0 || fabs (((double) grid5_xyz.z - load_min.xyz.z) * 0.533333) > 9.9)
        {
          *count = 1;
          load_point[0] = load_min;
        }


      /*  If the difference between the average depth and the maximum depth is greater than 9.9 fathoms (4800fps) we save the maximum location and depth.  */

      if (fabs (((double) grid5_xyz.z - load_max.xyz.z) * 0.533333) > 9.9)
        {
          load_point[*count] = load_min;
          (*count)++;
        }


      /*  If neither of the above cases was true we want to use the average depth.  */

      if (*count == 0)
        {
          *count = 1;
          load_point[0].coord.x = -1;
          load_point[0].xyz = grid5_xyz;
        }
    }


  /*  If there were no real or hand drawn depths we want to use the average depth.  */

  else
    {
      *count = 1;
      load_point[0].coord.x = -1;
      load_point[0].xyz = grid5_xyz;
    }
}


int32_t main (int32_t argc, char *argv[])
{
  char               c;
  extern char        *optarg;
  extern int         optind;
  int32_t            i, j, k, option_index = 0, chrtr2_handle[2], grid5_cols, grid5_rows, grid5_width, grid5_height;
  int32_t            grid_cols, grid_rows, load_count = 0, total_count = 0, row_filter, col_filter, percent = 0 , old_percent = -1;
  char               input_file[512], output_file[512];
  CHRTR2_HEADER      chrtr2_header[2];
  CHRTR2_RECORD      chrtr2_record[2];
  float              **grid5, *array, min_z, max_z;
  double             grid5_size_x, grid5_size_y;
  NV_F64_XYMBR       mbr5, new_mbr, mbr1;
  NV_F64_COORD3      xyz, grid5_xyz;
  NV_F64_COORD2      xy;
  LOAD_POINT         load_point[2], *total_load = NULL;
  NV_I32_COORD2      coord;


  printf ("\n\n %s \n\n\n", VERSION);


  /*  No command line options at present but you never know...  */

  while (NVTrue) 
    {
      static struct option long_options[] = {{0, no_argument, 0, 0}};

      c = (char) getopt_long (argc, argv, "h", long_options, &option_index);
      if (c == -1) break;

      switch (c) 
        {
        case 0:

          switch (option_index)
            {
            case 0:
              break;
            }
          break;


          /*  Placeholder.  */

        case 'h':
          break;

        default:
          usage ();
          break;
        }
    }


  /* Make sure we got the mandatory file name.  */
  
  if (optind >= argc) usage ();


  strcpy (input_file, argv[optind]);


  /*  Make the output file name.  */

  strcpy (output_file, input_file);
  sprintf (&output_file[strlen (output_file) - 4], "__filter5.ch2");

  fprintf (stderr, "Input file  : %s\n", input_file);
  fprintf (stderr, "Output file : %s\n\n", output_file);

 
  /*  Open the input file.  */

  chrtr2_handle[0] = chrtr2_open_file (input_file, &chrtr2_header[0], CHRTR2_READONLY);

  if (chrtr2_handle[0] < 0)
    {
      fprintf (stderr, "\n\nThe file %s is not a CHRTR2 file or there was an error reading the file.\nThe error message returned was:%s\n\n",
               input_file, chrtr2_strerror ());
      exit (-1);
    }


  /*  Define the MBR for the low resolution grid.  */

  grid5_size_x = chrtr2_header[0].lon_grid_size_degrees * 5.0;
  grid5_size_y = chrtr2_header[0].lat_grid_size_degrees * 5.0;
  grid5_width = NINT (chrtr2_header[0].width / 5.0);
  grid5_height = NINT (chrtr2_header[0].height / 5.0);

  mbr5.min_x = chrtr2_header[0].mbr.wlon;
  mbr5.min_y = chrtr2_header[0].mbr.slat;
  mbr5.max_x = mbr5.min_x + (double) grid5_width * grid5_size_x;
  mbr5.max_y = mbr5.min_y + (double) grid5_height * grid5_size_y;


  /*  Add the filter border to the MBR.  */

  mbr5.min_x -= ((double) FILTER * grid5_size_x);
  mbr5.min_y -= ((double) FILTER * grid5_size_y);
  mbr5.max_x += ((double) FILTER * grid5_size_x);
  mbr5.max_y += ((double) FILTER * grid5_size_y);


  /*  Number of rows and columns in the area  */

  grid5_rows = NINT ((mbr5.max_y - mbr5.min_y) / grid5_size_y);
  grid5_cols = NINT ((mbr5.max_x - mbr5.min_x) / grid5_size_x);


  /*  We're going to let MISP/SURF handle everything in zero based units of the bin size.  That is, we subtract off the
      west lon from longitudes then divide by the grid size in the X direction.  We do the same with the latitude using
      the south latitude.  This will give us values that range from 0.0 to grid5_cols in longitude and 0.0 to
      grid5_rows in latitude.  */

  new_mbr.min_x = 0.0;
  new_mbr.min_y = 0.0;
  new_mbr.max_x = (double) grid5_cols;
  new_mbr.max_y = (double) grid5_rows;


  misp_init (1.0, 1.0, 0.05, 4, 20.0, 20, 999999.0, -999999.0, -2, new_mbr);


  for (i = 0 ; i < chrtr2_header[0].height ; i++)
    {
      coord.y = i;

      for (j = 0 ; j < chrtr2_header[0].width ; j++)
        {
          coord.x = j;

          if (chrtr2_read_record (chrtr2_handle[0], coord, &chrtr2_record[0]))
            {
              chrtr2_perror ();
              exit (-1);
            }

          chrtr2_get_lat_lon (chrtr2_handle[0], &xy.y, &xy.x, coord);


          /*
              Load the points.

              IMPORTANT NOTE:  MISP and GMT (by default) grid using corner posts.  That is, the data in a bin is assigned to the 
              lower left corner of the bin.  Normal gridding/binning systems use the center of the bin.  Because of this we need
              to lie to MISP/GMT and tell them that the point is really half a bin lower and to the left.  This is extremely
              confusing but it works ;-)
          */


          xyz.x = (xy.x - mbr5.min_x) / grid5_size_x;
          xyz.y = (xy.y - mbr5.min_y) / grid5_size_y;
          xyz.z = chrtr2_record[0].z;

          misp_load (xyz);
        }

      percent = NINT (((float) i / (float) chrtr2_header[0].height) * 100.0);
      if (percent != old_percent)
        {
          fprintf (stderr, "Loading data for low resolution grid - %03d%% complete\r", percent);
          old_percent = percent;
        }
    }


  fprintf (stderr, "                                                                   \r");
  fprintf (stderr, "\nLow resolution grid load complete\n\n");


  fprintf (stderr, "Processing low resolution grid\n");


  misp_proc ();


  fprintf (stderr, "Processing low resolution grid complete\n");


  /*  Allocating one more than grid5_cols and grid5_rows due to constraints of old chrtr (see comments in misp_funcs.c)  */

  grid5 = (float **) malloc ((grid5_rows + 1) * sizeof (float *));
  if (grid5 == NULL)
    {
      perror ("Allocating grid5 in main.c");
      exit (-1);
    }

  for (i = 0 ; i <= grid5_rows ; i++)
    {
      grid5[i] = (float *) malloc ((grid5_cols + 1) * sizeof (float));
      if (grid5[i] == NULL)
        {
          perror ("Allocating grid5[i] in main.c");
          exit (-1);
        }
    }


  for (i = 0 ; i <= grid5_rows ; i++)
    {
      if (!misp_rtrv (grid5[i])) break;


      percent = NINT (((float) i / (float) grid5_rows) * 100.0);
      if (percent != old_percent)
        {
          fprintf (stderr, "Retrieving data for low resolution grid - %03d%% complete\r", percent);
          old_percent = percent;
        }
    }

  fprintf (stderr, "                                                                   \r");
  fprintf (stderr, "\nLow resolution grid retrieval complete\n\n");



  /*  Now retrieve the points to be loaded by comparing the low resolution grid to the original file.  */

  for (i = FILTER ; i < grid5_rows ; i++)
    {
      grid5_xyz.y = mbr5.min_y + (((double) i) * grid5_size_y);

      for (j = FILTER ; j < grid5_cols ; j++)
        {
          grid5_xyz.x = mbr5.min_x + (((double) j) * grid5_size_x);

          grid5_xyz.z = grid5[i][j];

          get_load_values (grid5_xyz, i, j, chrtr2_handle[0], chrtr2_header[0], load_point, &load_count);

          for (k = 0 ; k < load_count ; k++)
            {
              total_load = (LOAD_POINT *) realloc (total_load, (total_count + 1) * sizeof (LOAD_POINT));

              if (total_load == NULL)
                {
                  perror ("Allocating total_load in main.c");
                  exit (-1);
                }

              total_load[total_count] = load_point[k];
              total_count++;
            }
        }

      percent = NINT (((float) i / (float) grid5_rows) * 100.0);
      if (percent != old_percent)
        {
          fprintf (stderr, "Comparing low resolution grid to original grid - %03d%% complete\r", percent);
          old_percent = percent;
        }
    }


  fprintf (stderr, "                                                                                   \r");
  fprintf (stderr, "\nLow resolution grid to original grid comparison complete\n\n");


  /*  At this point we can safely free the low resolution grid since we've got all the information we need from it.  */

  for (i = 0 ; i <= grid5_rows ; i++) free (grid5[i]);
  free (grid5);


  /*  Now we get to make a new MBR for the new file at the original resolution.  */

  mbr1.min_x = chrtr2_header[0].mbr.wlon;
  mbr1.min_y = chrtr2_header[0].mbr.slat;
  mbr1.max_x = chrtr2_header[0].mbr.elon;
  mbr1.max_y = chrtr2_header[0].mbr.nlat;


  /*  Add the filter border to the MBR  */

  mbr1.min_x -= ((double) FILTER * chrtr2_header[0].lon_grid_size_degrees);
  mbr1.min_y -= ((double) FILTER * chrtr2_header[0].lat_grid_size_degrees);
  mbr1.max_x += ((double) FILTER * chrtr2_header[0].lon_grid_size_degrees);
  mbr1.max_y += ((double) FILTER * chrtr2_header[0].lat_grid_size_degrees);


  /*  Number of rows and columns in the area  */

  grid_cols = NINT ((mbr1.max_x - mbr1.min_x) / chrtr2_header[0].lon_grid_size_degrees);
  grid_rows = NINT ((mbr1.max_y - mbr1.min_y) / chrtr2_header[0].lat_grid_size_degrees);


  row_filter = grid_rows - FILTER;
  col_filter = grid_cols - FILTER;


  /*  And, once again, make a phony new MBR with units of 1.0  */

  new_mbr.min_x = 0.0;
  new_mbr.min_y = 0.0;
  new_mbr.max_x = (double) grid_cols;
  new_mbr.max_y = (double) grid_rows;


  misp_init (1.0, 1.0, 0.05, 4, 20.0, 20, 999999.0, -999999.0, -2, new_mbr);


  for (i = 0 ; i < total_count ; i++)
    {
      /*
        Load the points.

        IMPORTANT NOTE:  MISP and GMT (by default) grid using corner posts.  That is, the data in a bin is assigned to the 
        lower left corner of the bin.  Normal gridding/binning systems use the center of the bin.  Because of this we need
        to lie to MISP/GMT and tell them that the point is really half a bin lower and to the left.  This is extremely
        confusing but it works ;-)
      */

      xyz.x = (total_load[i].xyz.x - mbr1.min_x) / chrtr2_header[0].lon_grid_size_degrees;
      xyz.y = (total_load[i].xyz.y - mbr1.min_y) / chrtr2_header[0].lat_grid_size_degrees;
      xyz.z = total_load[i].xyz.z;

      misp_load (xyz);


      percent = NINT (((float) i / (float) total_count) * 100.0);
      if (percent != old_percent)
        {
          fprintf (stderr, "Loading data for final grid - %03d%% complete\r", percent);
          old_percent = percent;
        }
    }


  fprintf (stderr, "                                                                   \r");
  fprintf (stderr, "\nFinal grid load complete\n\n");


  fprintf (stderr, "Processing final grid\n");


  misp_proc ();


  fprintf (stderr, "Processing final grid complete\n");


  /*  Allocating one more than grid_cols due to constraints of old chrtr (see comments in misp_funcs.c)  */

  array = (float *) malloc ((grid_cols + 1) * sizeof (float));

  if (array == NULL)
    {
      perror ("Allocating array in main.c");
      exit (-1);
    }


  /*  Try to create and open the new chrtr2 file.  */

  chrtr2_header[1] = chrtr2_header[0];
  strcpy (chrtr2_header[1].creation_software, VERSION);
  chrtr2_handle[1] = chrtr2_create_file (output_file, &chrtr2_header[1]);
  if (chrtr2_handle[1] < 0)
    {
      chrtr2_perror ();
      exit (-1);
    }


  min_z = 9999999999.0;
  max_z = -9999999999.0;


  /*  This is where we stuff the new interpolated surface in to the new CHRTR2.  */

  for (i = 0 ; i < grid_rows ; i++)
    {
      if (!misp_rtrv (array)) break;


      /*  Only use data that aren't in the filter border  */

      if (i >= FILTER && i < row_filter)
        {
          coord.y = i - FILTER;

          for (j = 0 ; j < grid_cols ; j++)
            {
              /*  Only use data that aren't in the filter border  */

              if (j >= FILTER && j < col_filter)
                {
                  coord.x = j - FILTER;


                  /*  Make sure we're inside the CHRTR2 bounds.  */

                  if (coord.y >= 0 && coord.y < chrtr2_header[1].height && coord.x >= 0 && coord.x < chrtr2_header[1].width)
                    {
                      /*  Read the record so that we can preserve all of the rest of the information from the original file
                          (i.e. status, uncertainty...)  */

                      chrtr2_read_record (chrtr2_handle[0], coord, &chrtr2_record[0]);


                      /*  Stuff everything from the original record into the new record with the exception of the Z value.  */

                      chrtr2_record[1] = chrtr2_record[0];
                      chrtr2_record[1].z = array[j];


                      min_z = MIN (chrtr2_record[1].z, min_z);
                      max_z = MAX (chrtr2_record[1].z, max_z);

                      chrtr2_write_record (chrtr2_handle[1], coord, chrtr2_record[1]);
                    }
                }
            }
        }

      percent = NINT (((float) i / (float) grid_rows) * 100.0);
      if (percent != old_percent)
        {
          fprintf (stderr, "Retrieving data for final grid - %03d%% complete\r", percent);
          old_percent = percent;
        }
    }


  fprintf (stderr, "                                                                   \r");
  fprintf (stderr, "\nFinal grid retrieval complete\n\n");


  /*  Now we want to read through the high and low points and replace the interpolated values in the new file
      with the values from the original file.  */

  for (i = 0 ; i < total_count ; i++)
    {
      /*  If the coordinates aren't set properly, this was an interpolated point anyway.  */

      if (total_load[i].coord.x >= 0)
        {
          chrtr2_read_record (chrtr2_handle[0], total_load[i].coord, &chrtr2_record[0]);
          chrtr2_read_record (chrtr2_handle[1], total_load[i].coord, &chrtr2_record[1]);


          chrtr2_record[1].z = chrtr2_record[0].z;

          min_z = MIN (chrtr2_record[1].z, min_z);
          max_z = MAX (chrtr2_record[1].z, max_z);

          chrtr2_write_record (chrtr2_handle[1], total_load[i].coord, chrtr2_record[1]);
        }
    }


  free (total_load);


  chrtr2_close_file (chrtr2_handle[0]);


  /*  Update the header with the observed min and max values.  */

  chrtr2_header[1].min_observed_z = min_z;
  chrtr2_header[1].max_observed_z = max_z;

  chrtr2_update_header (chrtr2_handle[1], chrtr2_header[1]);

  chrtr2_close_file (chrtr2_handle[1]);


  fprintf (stderr, "\n\n%s complete\n\n\n", argv[0]);


  /*  Please ignore the following line.  It is useless.  Except...

      On some versions of Ubuntu, if I compile a program that doesn't use the math
      library but it calls a shared library that does use the math library I get undefined
      references to some of the math library functions even though I have -lm as the last
      library listed on the link line.  This happens whether I use qmake to build the
      Makefile or I have a pre-built Makefile.  Including math.h doesn't fix it either.
      The following line forces the linker to bring in the math library.  If there is a
      better solution please let me know at area.based.editor AT gmail DOT com.  */

  float ubuntu; ubuntu = 4.50 ; ubuntu = fmod (ubuntu, 1.0);


  return (0);
}
