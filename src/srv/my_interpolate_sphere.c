/*
   interpolate_sphere.c

   To sphere-to-sphere interpolation to interpolate a surface object
   at the vertices of a sphere in standardized space.

   interpolate_sphere surf.obj inflate.obj sphere.obj out.obj 

   Values: surf.obj = input object file
           inflate.obj = input object file inflated to a sphere
           sphere.obj = stereotaxic sphere model
           out.obj = output object file

   By: Claude Lepage, May 2011.

   COPYRIGHT: McConnell Brain Imaging Cenbiter, 
              Montreal Neurological Institute,
              Department of Psychology,
              McGill University, Montreal, Quebec, Canada. 
  
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
*/

#include <math.h>
#include <stdio.h>
#include "minc2.h"
#include "volume_io.h"
#include "bicpl.h"


// Prototypes of functions in this file.

static void usage( char * );
static VIO_Status read_surface_obj( char*, int *, VIO_Point *[],
                                VIO_Vector *[], int *, int *[], int **, int *** );
static VIO_Status get_surface_neighbours( polygons_struct *, int **,
                                      int ***  );
static void save_surface_obj( char*, int, VIO_Point *, VIO_Vector *, int, int []);

static void compute_triangle_angles( int, int *, VIO_Point * );
static void compute_edge_stats( int, int *, VIO_Point * );
static void compute_triangle_normal( VIO_Point, VIO_Point, VIO_Point, VIO_Real[3] );
void smooth( int, int, int, VIO_Real, int *, VIO_Point *, float );

// Main program.
/*int main(int argc, char** argv)
{

    int* n_ngh;
    int** ngh;
    interpolate_sphere( argv[1], &n_ngh, &ngh);


return 0;
}*/


int interpolate_sphere( char * input_object_filename, int** n_ngh, int*** ngh ) {

  int      i, j, k, jj, kk, pp, v1, v2, opp1, opp2, t1, t2, found3;
  int      target_nodes, changed, num_swapped, histo[101];
  VIO_Real     thresholdlen, minlen, maxlen, factor;

  int      n_points, n_points_1;  // number of grid points per object
  int      n_elems, n_elems_1;    // number of triangles per object
  VIO_Point  * original_coords;       // coordinates
  VIO_Point  * coords;                // coordinates
  VIO_Vector * normals;               // normal vectors
  int    * connec;                // connectivity

 

  // Read the surface.
  if( read_surface_obj( input_object_filename, &n_points, &coords, &normals,
                        &n_elems, &connec, n_ngh, ngh ) != VIO_OK ) {
      
    return 1;
  }


  return 0;
}




// -------------------------------------------------------------------
// Help message on how to use this module.
//
static void usage( char * executable_name ) {

  char*  usage_format = "\
Usage: %s surf.obj inflate.obj sphere.obj out.obj\n\
Values: surf.obj = input object file\n\
        inflate.obj = input object file inflated to a sphere\n\
        sphere.obj = stereotaxic sphere model\n\
        out.obj = output object file\n\n";

  print_error( usage_format, executable_name );
}


// -------------------------------------------------------------------
// Load the cortical surface.
//
// filename: name of the .obj file
// n_points: the number of the vertices
// points: (x,y,z) coordinates
// normals: normal vectors
// n_elem: number of triangles
// connec: connectivity of triangles
// n_neighbours: number of vertices around each node
// neighbours: the set of ordered triangle consisting of the vertices
//
static VIO_Status read_surface_obj( char* filename,
                                 int * n_points,
                                 VIO_Point * points[],
                                 VIO_Vector * normals[],
                                 int * n_elem,
                                 int * connec[],
                                 int ** n_neighbours,
                                 int *** neighbours ) {

  int               i, n_objects;
  object_struct  ** object_list;
  polygons_struct * surface;
  VIO_File_formats*      format;
  char*            expanded;

  expanded = expand_filename( filename );   // why?????

  int err = input_graphics_file( expanded, &format, &n_objects,
                                 &object_list );

  if( err != VIO_OK ) {
    print_error( "Error reading file %s\n", expanded );
    return( VIO_ERROR );
  }

  if( /*n_objects !=   || commented outbecause i dont know what this does*/
      ( n_objects == 1 && get_object_type(object_list[0]) != POLYGONS ) ) {
    print_error( "Error in contents of file %s\n", expanded );
    return( VIO_ERROR );
  }

  delete_string( expanded );

  surface = get_polygons_ptr( object_list[0] );

  int ntri = 0, nquad = 0, unknown = 0;
  int start_ind = 0;
  for( i = 0; i < surface->n_items; i++ ) {
    int nn = surface->end_indices[i] - start_ind;
    start_ind = surface->end_indices[i];
    if( nn == 3 ) {
      ntri++;
    } else {
     if( nn == 4 ) {
       nquad++;
     } else {
       unknown++;
       //printf( "face with %d nodes\n", nn );
     }
   }
  }
  //printf( "%d triangles, %d quads, %d unknown faces in mesh\nsurface points %d\n", ntri, nquad, unknown, surface->n_points );

  // Check if all polygons are triangles.

  if( 3 * surface->n_items != surface->end_indices[surface->n_items-1] ) {
    printf( "Error: Surface must contain only triangular polygons.\n" );
    delete_object_list( n_objects, object_list );
    return VIO_ERROR;
  }

  // Make a copy of the coordinates, the normals, and the
  // connectivity since delete_object_list will destroy them.

  *n_points = surface->n_points;
  //printf("npoints %d\n", *n_points);
  *n_elem = surface->n_items;
  ALLOC( *points, surface->n_points );
  if( !(*points) ) {
    printf( "Error allocating memory for points.\n" );
    exit( 1 );
  }
  ALLOC( *normals, surface->n_points );
  if( !(*normals) ) {
    printf( "Error allocating memory for normals.\n" );
    exit( 1 );
  }

  for( i = 0; i < *n_points; i++ ) {
    (*points)[i].coords[0] = surface->points[i].coords[0];
    (*points)[i].coords[1] = surface->points[i].coords[1];
    (*points)[i].coords[2] = surface->points[i].coords[2];
    (*normals)[i].coords[0] = surface->normals[i].coords[0];
    (*normals)[i].coords[1] = surface->normals[i].coords[1];
    (*normals)[i].coords[2] = surface->normals[i].coords[2];
  }

  if( connec ) {
    ALLOC( *connec, 3*surface->n_items );
    if( !(*connec) ) {
      printf( "Error allocating memory for connec.\n" );
      exit( 1 );
    }
    for( i = 0; i < 3*surface->n_items; i++ ) {
      (*connec)[i] = surface->indices[i];
    }
  }
//printf("okay so far %d %d\n", *n_neighbours, neighbours);
  //if( n_neighbours && neighbours ) {
    //printf("got to here\n");
      get_surface_neighbours( surface, n_neighbours, neighbours );
 // }

  delete_object_list( n_objects, object_list );

  return( VIO_OK );
}


// -------------------------------------------------------------------
// Construct the edges around each node. The edges are sorted to
// make an ordered closed loop.
//
static VIO_Status get_surface_neighbours( polygons_struct * surface,
                                       int ** n_ngh,
                                       int *** ngh ) {

  int    i, j, k, jj;
  int  * tri;
  int  * ngh_array;

  // Check if all polygons are triangles.

  if( 3 * surface->n_items != surface->end_indices[surface->n_items-1] ) {
    printf( "Surface must contain only triangular polygons.\n" );
    return VIO_ERROR;
  }

  // Check if the node numbering starts at 0 or 1.

  int min_idx, max_idx;
  min_idx = 100*surface->n_points;  // anything big
  max_idx = 0;                      // anything small

  for( i = 0; i < 3*surface->n_items; i++ ) {
    if( surface->indices[i] < min_idx ) min_idx = surface->indices[i];
    if( surface->indices[i] > max_idx ) max_idx = surface->indices[i];
  }

  // Shift numbering to start at zero, for array indexing. Note
  // that we don't care if surface->indices array is modified.

  if( min_idx != 0 ) {
    for( i = 0; i < 3*surface->n_items; i++ ) {
      surface->indices[i] -= min_idx;
    }
  }

  // Count number of triangles attached to each node.

  ALLOC( *n_ngh, surface->n_points );
  if( !(*n_ngh) ) {
    printf( "Error allocating memory for n_ngh.\n" );
    exit( 1 );
  }
  ALLOC( *ngh, surface->n_points );
  if( !(*ngh) ) {
    printf( "Error allocating memory for ngh.\n" );
    exit( 1 );
  }
  ALLOC( ngh_array, 3*surface->n_items );
  if( !ngh_array ) {
    printf( "Error allocating memory for ngh_array.\n" );
    exit( 1 );
  }

  for( i = 0; i < surface->n_points; i++ ) {
    (*n_ngh)[i] = 0;
  }

  for( i = 0; i < 3*surface->n_items; i++ ) {
    (*n_ngh)[surface->indices[i]]++;
    ngh_array[i] = -1;
  }

  int max_ngh = 0;
  int sum_ngh = 0;
  for( i = 0; i < surface->n_points; i++ ) {
    (*ngh)[i] = &(ngh_array[sum_ngh]);
    sum_ngh += (*n_ngh)[i];
    max_ngh = MAX( max_ngh, (*n_ngh)[i] );
  }

  // At first, store the indices of the triangles in the neighbours.
  for( i = 0; i < surface->n_items; i++ ) {
    for( j = 0; j < 3; j++ ) {
      jj = surface->indices[3*i+j];
      for( k = 0; k < (*n_ngh)[jj]; k++ ) {
        if( (*ngh)[jj][k] == -1 ) {
          (*ngh)[jj][k] = i;
          break;
        }
      }
    }
  }

  // Now create a sort closed loop of the node neighbours.
  // This is needed by the parametric=0 FEM algorithm.
  //
  //         1 ----- 2
  //          /\   /\
  //         /  \ /  \
  //       0 ----P---- 3
  //         \  / \  /
  //          \/   \/
  //         5 ----- 4
  //

  int * tmp;
  ALLOC( tmp, 2*max_ngh );
  if( !tmp ) {
    printf( "Error allocating memory for tmp.\n" );
    exit( 1 );
  }


   for( i = 0; i < surface->n_points; i++ ) {
    for( k = 0; k < (*n_ngh)[i]; k++ ) {
      tri = &(surface->indices[3*(*ngh)[i][k]]);
      for( j = 0; j < 3; j++ ) {
        if( tri[j] == i ) break;
      }
      tmp[2*k+0] = tri[(j+1)%3];
      tmp[2*k+1] = tri[(j+2)%3];
    }

    (*ngh)[i][0] = tmp[0];
    (*ngh)[i][1] = tmp[1];
    for( k = 2; k < (*n_ngh)[i]; k++ ) {
      for( j = 1; j < (*n_ngh)[i]; j++ ) {
        if( tmp[2*j] == (*ngh)[i][k-1] || tmp[2*j+1] == (*ngh)[i][k-1] ) {
          if( tmp[2*j] == (*ngh)[i][k-1] ) {
            (*ngh)[i][k] = tmp[2*j+1];
          } else {
            (*ngh)[i][k] = tmp[2*j];
          }
          tmp[2*j] = -1;
          tmp[2*j+1] = -1;
          break;
        }
      }
    }
  }

//  n_neighbours_return = *n_ngh;
//  neighbours_return = *ngh;
 /*for( i = 0; i < surface->n_points; i++ ) 
  {
    printf("%d\n");
    for( k = 0; k < (*n_ngh)[i] ; k++ )
        printf("\t%d\n", (*ngh)[i][k]  );
  }*/


  FREE( tmp );

  return VIO_OK;

}


