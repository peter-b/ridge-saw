/*
 * Surrey Space Centre ridge/self-avoiding walk tool
 * Copyright (C) 2012  Peter Brett <p.brett@surrey.ac.uk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"

#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <glib.h>
#include <gsl/gsl_randist.h>

#define GETOPT_OPTIONS "i:r::d:t:n:s:h"

#include <ridgeutil.h>
#include <ridgeio.h>

enum GenerateMode {
  GENERATE_SPECKLE = 0,
  GENERATE_NORM,
};

void
usage (char *name, int status)
{
  printf (
"Usage: %s OPTION... [FILE]\n"
"\n"
"Generate ridge data for self-avoiding walk analysis.\n"
"\n"
"  -r [TYPE]       Generate random image data [default: S]\n"
"  -d SIZE         Size for random tiles [default: 2048]\n"
"  -t SCALE        Ridge detection scale [default: 0]\n"
"  -n NUM          Target data point count for random generation\n"
"  -s SEED         Random seed.\n"
"  -h              Display this message and exit\n"
"\n"
"Detect ridge lines and output step count and end-to-end distance for\n"
"comparison with self-avoiding walk statistics.  Two modes are\n"
"available:\n"
"\n"
"  - If a FILE was specified, image data is loaded from FILE, and\n"
"    the number of data points is determined automatically.\n"
"\n"
"  - If the '-r' option was given, random noise images are generated\n"
"    and used to obtain line data.  The '-r' option controls the\n"
"    noise function used; the TYPE must be 'S' (default) or 'N'.  The\n"
"    '-d' option controls how large the generated images are.  If the\n"
"    '-n' option is given, images will be repeatedly generated until\n"
"    NUM data points have been created.  The '-s' option allows the\n"
"    random number generator seed to be overridden.\n"
"\n"
"If an OUTFILE was specified, CSV data is output to that file;\n"
"otherwise, output is to standard output.\n"
"\n"
"The RIDGETOOL environment variable can be set to control the path to\n"
"the 'ridgetool' program.\n"
"\n"
"Please report bugs to %s.\n",
name, PACKAGE_BUGREPORT);
  exit (status);
}

RioData *
run_ridgetool_get_data (const char *filename, float scale)
{
  g_assert (filename);

  /* Figure out how to call ridgetool */
  const gchar *ridgetool_path = g_getenv ("RIDGETOOL");
  if (ridgetool_path == NULL) ridgetool_path = "ridgetool";

  /* Get a temporary filename. FIXME we don't use this in a secure
   * way, unfortunately. */
  gchar *tmpfile = g_strdup ("ridge-saw.XXXXXX");
  int tmpfd = mkstemp (tmpfile);
  if (tmpfd == -1) {
    const char *msg = errno ? strerror (errno) : "Unexpected error";
    fprintf (stderr, "ERROR: Failed to create temporary file: %s\n\n",
             msg);
    exit (2);
  }

  /* Format scale as string */
  gchar *sscale = g_strdup_printf ("-t%f", scale);

  /* Invoke ridgetool */
  const gchar *argv[] = {ridgetool_path,
                         "-l",
                         sscale,
                         "-i0",
                         filename,
                         tmpfile,
                         NULL};
  gint exit_status;
  gchar *err_output = NULL;
  GError *err = NULL;
  g_spawn_sync (NULL /* working dir */,
                (gchar **) argv,
                NULL /* envp */,
                G_SPAWN_SEARCH_PATH /* flags */,
                NULL /* child setup */,
                NULL /* child setup data */,
                NULL /* standard output */,
                &err_output /* standard error */,
                &exit_status,
                &err);
  if (err != NULL) {
    fprintf (stderr, "ERROR: Failed to run '%s': %s\n\n",
             ridgetool_path, err->message);
    exit (3);
  }
  if (exit_status != 0) {
    fprintf (stderr, "ERROR: '%s' failed:\n%s\n\n",
             ridgetool_path, err_output);
    exit(3);
  }
  g_free (sscale);

  /* Load ridge data */
  RioData *data = rio_data_from_file (tmpfile);
  if (data == NULL) {
    const char *msg = errno ? strerror (errno) : "Unexpected error";
    fprintf (stderr, "ERROR: Failed to load ridge data from '%s': %s\n\n",
             tmpfile, msg);
    exit (2);
  }

  /* Clean up */
  close (tmpfd);
  unlink (tmpfile);
  g_free (tmpfile);

  return data;
}

int
dump_saw_stats (RioData *data, FILE *fp)
{
  g_assert (data);
  g_assert (fp);
  g_assert (rio_data_get_type (data) == RIO_DATA_LINES);

  for (int i = 0; i < rio_data_get_num_entries (data); i++) {
    RioLine *l = rio_data_get_line (data, i);
    int len = rio_line_get_length (l);
    RioPoint *start = rio_line_get_point (l, 0);
    RioPoint *end = rio_line_get_point (l, len - 1);

    /* Calculate distance */
    double start_row, start_col, end_row, end_col;
    rio_point_get_subpixel (start, &start_row, &start_col);
    rio_point_get_subpixel (end, &end_row, &end_col);
    double dx = floor (end_col) - floor (start_col);
    double dy = floor (end_row) - floor (start_row);
    double dist = sqrt (dx*dx + dy*dy);

    /* Output is in the format "num_steps, distance" */
    if (fprintf (fp, "%i, %f\n", len - 1, dist) < 0) {
      return 0;
    }
  }
  return 1; /* Success */
}


int
main (int argc, char **argv)
{
  int gen_mode = -1;
  int gen_size = 2048;
  int gen_target = -1;
  int gen_seed = -1;
  char *infile = NULL;
  char *outfile = NULL;
  float scale = 0;
  int c, status;

  /* Parse command-line arguments */
  while ((c = getopt (argc, argv, GETOPT_OPTIONS)) != -1) {
    switch (c) {
    case 'i':
      if (gen_mode != -1) {
        fprintf (stderr,
                 "ERROR: Only one of '-i' or '-r' options may be given.\n\n");
        usage (argv[0], 1);
      }
      infile = optarg;
      break;
    case 'r':
      if (infile != NULL) {
        fprintf (stderr,
                 "ERROR: Only one of '-i' or '-r' options may be given.\n\n");
        usage (argv[0], 1);
      }
      if (optarg == NULL) {
        gen_mode = GENERATE_SPECKLE;
      } else {
        switch (optarg[0]) {
        case 'S': gen_mode = GENERATE_SPECKLE;
        case 'N': gen_mode = GENERATE_NORM;
        default:
          fprintf (stderr, "ERROR: Bad argument '%s' to -r option.\n\n",
                   optarg);
          usage (argv[0], 1);
        }
      }
      break;
    case 'd':
      status = sscanf (optarg, "%i", &gen_mode);
      if (status != 1 || gen_mode < 1) {
        fprintf (stderr, "ERROR: Bad argument '%s' to -d option.\n\n",
                 optarg);
        usage (argv[0], 1);
      }
      break;
    case 't':
      status = sscanf (optarg, "%f", &scale);
      if (status != 1 || scale < 0) {
        fprintf (stderr, "ERROR: Bad argument '%s' to -t option.\n\n",
                 optarg);
        usage (argv[0], 1);
      }
      break;
    case 'n':
      status = sscanf (optarg, "%i", &gen_target);
      if (status != 1 || gen_target < 1) {
        fprintf (stderr, "ERROR: Bad argument '%s' to -n option.\n\n",
                 optarg);
        usage (argv[0], 1);
      }
      break;
    case 's':
      status = sscanf (optarg, "%i", &gen_seed);
      if (status != 1 || gen_seed < 1) {
        fprintf (stderr, "ERROR: Bad argument '%s' to -s option.\n\n",
                 optarg);
        usage (argv[0], 1);
      }
      break;
    case 'h':
      usage (argv[0], 0);
      break;
    case '?':
      if ((optopt != ':') && (strchr (GETOPT_OPTIONS, optopt) != NULL)) {
        fprintf (stderr, "ERROR: -%c option requires an argument.\n\n", optopt);
      } else if (isprint (optopt)) {
        fprintf (stderr, "ERROR: Unknown option -%c.\n\n", optopt);
      } else {
        fprintf (stderr, "ERROR: Unknown option character '\\x%x'.\n\n",
                 optopt);
      }
      usage (argv[0], 1);
    default:
      g_assert_not_reached ();
    }
  }

  if (gen_mode == -1 && !outfile) {
    fprintf (stderr, "ERROR: You must specify '-r' or '-i' options.\n\n");
    usage (argv[0], 1);
  }

  FILE *outfp = stdout;
  if (outfile != NULL) {
    outfp = fopen (outfile, "wb");
    if (outfp == NULL) {
      const char *msg = errno ? strerror (errno) : "Unexpected error";
      fprintf (stderr, "ERROR: Failed to open output file '%s': %s\n\n",
               outfile, msg);
      exit (4);
    }
  }

  if (infile != NULL) {
    /* Load and process input file */
    RioData *data = run_ridgetool_get_data (infile, scale);
    status = dump_saw_stats (data, outfp);
    if (!status) {
      const char *msg = errno ? strerror (errno) : "Unexpected error";
      fprintf (stderr, "ERROR: Output failed: %s\n\n", msg);
      exit (4);
    }
    rio_data_destroy (data);

  } else if (gen_mode != -1) {
    int N = 0;

    /* Initialise RNG. */
    gsl_rng_env_setup ();
    gsl_rng *rng = gsl_rng_alloc (gsl_rng_default);
    if (gen_seed >= 0) {
      gsl_rng_set (rng, (unsigned int) gen_seed);
    }
    fprintf (stderr, "Random number seed: %lu (%s)\n",
             (gen_seed >= 0) ? gen_seed : gsl_rng_default_seed,
             gsl_rng_name (rng));

    /* Get a temporary filename. FIXME we don't use this in a secure
     * way, unfortunately. */
    gchar *tmpfile = g_strdup ("ridge-saw.XXXXXX");
    int tmpfd = mkstemp (tmpfile);
    if (tmpfd == -1) {
      const char *msg = errno ? strerror (errno) : "Unexpected error";
      fprintf (stderr, "ERROR: Failed to create temporary file: %s\n\n",
               msg);
      exit (5);
    }

    /* Repeatedly generate and process random images */
    RutSurface *img = rut_surface_new (gen_size, gen_size);
    do {

      /* Generate random data */
      for (int i = 0; i < img->rows; i++) {
        for (int j = 0; j < img->cols; j++) {
          double val;
          switch (gen_mode) {
          case GENERATE_NORM:
            val = gsl_ran_gaussian (rng, 1);
            break;
          case GENERATE_SPECKLE:
            val = gsl_ran_rayleigh (rng, 1);
            break;
          default:
            g_assert_not_reached ();
          }
          RUT_SURFACE_REF(img,i,j) = (float) val;
        }
      }

      /* Output to TIFF file */
      status = rut_surface_to_tiff (img, tmpfile);
      if (!status) {
        fprintf (stderr, "ERROR: Failed to write image data to '%s'.\n\n",
                 tmpfile);
        exit (5);
      }

      /* Process TIFF file */
      RioData *data = run_ridgetool_get_data (tmpfile, scale);
      status = dump_saw_stats (data, outfp);
      if (!status) {
        const char *msg = errno ? strerror (errno) : "Unexpected error";
        fprintf (stderr, "ERROR: Output failed: %s\n\n", msg);
        exit (4);
      }
      rio_data_destroy (data);

    } while (N < gen_target);

    close (tmpfd);
    unlink (tmpfile);
    g_free (tmpfile);
    rut_surface_destroy (img);
    gsl_rng_free (rng);

  } else {
    g_assert_not_reached ();
  }

  if (outfp != stdout) {
    status = fclose (outfp);
    if (status != 0) {
      const char *msg = errno ? strerror (errno) : "Unexpected error";
      fprintf (stderr, "ERROR: Failed to close output file '%s': %s\n\n",
               outfile, msg);
      exit (4);
    }
  }
  exit (0);
}
