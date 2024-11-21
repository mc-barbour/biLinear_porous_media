/***************************************************************/

/* UDF for creating spatially varying porous media model coefficients

/* Created by Mike Barbour */
/* user must create a UDM for visualization*/

/* Now using a kd-tree search*/
/******************************************************/
#include "assert.h"
#include "udf.h"
#include "kdtree.h"
#include "malloc.h"
#include "bilinear_parameters.h"


DEFINE_PROFILE(inertial_res, vol_thread, i)
{
#if !RP_HOST
    real current_time;
    current_time = RP_Get_Real("flow-time");
 if (current_time == 0.0)   /* replace first_iteration */
 {

    real x[ND_ND];
    real x_f[ND_ND];
    real (*x_f_array)[ND_ND];
    real (*x_f_array_total)[ND_ND];

    cell_t c;
    face_t f;
    int k, jj;
    int tree_count = 0;

    double phi, perm, c2;
    double dist;
    double slope = (core_porosity - wall_porosity) / dist_threshold;

    int face_int, total_num_faces;

    void *kd, *kd_result, *kd_node;
    double pos[ND_ND];   
    
    kd = kd_create(ND_ND);

    Domain *domain = Get_Domain(aneurysm_vol_id);

    int face_ID = aneurysm_surf_id;
    Thread *boundary_thread = Lookup_Thread(domain, face_ID);
    int n_faces = THREAD_N_ELEMENTS_INT(boundary_thread);
    int n_faces2 = THREAD_N_ELEMENTS_INT(boundary_thread);
    

    /* Send number of faces to node 0 first time*/
    #if RP_NODE
        if(! I_AM_NODE_ZERO_P) /* Nodes 1,2... send the number of faces */
        {
            PRF_CSEND_INT(node_zero, &n_faces2, 1, myid);
        }
    #endif

    /* recieve number of faces at node 0 to compute total number of faces*/
    #if RP_NODE
        if(I_AM_NODE_ZERO_P)
        {
            total_num_faces = n_faces2;
            compute_node_loop_not_zero(k)
            {
                PRF_CRECV_INT(k, &n_faces2, 1, k);
                total_num_faces = total_num_faces + n_faces2;
            }
        }
    #endif


    
    /* Send number of face elements to node zero second time  */
    #if RP_NODE
        if(! I_AM_NODE_ZERO_P) /* Nodes 1,2... send the number of faces */
        {
            PRF_CSEND_INT(node_zero, &n_faces, 1, myid);
        }
    #endif

    /* allocate arrays */
    x_f_array = (real (*)[ND_ND])malloc(ND_ND*n_faces*sizeof(real));
    x_f_array_total = (real (*)[ND_ND])malloc(ND_ND*n_faces*sizeof(real));
    

    /* Create arrays of face coordinates at each node*/
    int noface = 0;
    begin_f_loop(f, boundary_thread)
    {
        F_CENTROID(x_f_array[f], f, boundary_thread);
        noface++;
    }
    end_f_loop(f, boundary_thread)

    /* Send face coordinates to node-zero*/
    #if RP_NODE
        if(! I_AM_NODE_ZERO_P) /* Only SEND data from nodes 1,2... */
        {
            PRF_CSEND_REAL(node_zero, x_f_array[0], ND_ND*n_faces, myid);
        }
    #endif

    /* Collect all face coordinate data at node zero */
    #if RP_NODE
        if(I_AM_NODE_ZERO_P)
        {
            x_f_array_total = (real (*)[ND_ND])malloc(ND_ND*total_num_faces*sizeof(real));
            
            /* Fill array with faces from node zero */
            begin_f_loop(f, boundary_thread)
            {            
                F_CENTROID(x_f_array_total[f], f, boundary_thread);
            }
            end_f_loop(f, boundary_thread)

            /* Fill array with faces from non-zero nodes */
            compute_node_loop_not_zero(k)
            {
    
                PRF_CRECV_INT(k, &n_faces, 1, k);

                /* Reallocate memory for arrays for node-i */               
                x_f_array=(real(*)[ND_ND])realloc(x_f_array,ND_ND*n_faces*sizeof(real));

                PRF_CRECV_REAL(k, x_f_array[0], ND_ND*n_faces, k);

                for (jj = 0; jj<n_faces; jj++){                    
                    x_f_array_total[noface+jj][0] = x_f_array[jj][0];
                    x_f_array_total[noface+jj][1] = x_f_array[jj][1];
                    x_f_array_total[noface+jj][2] = x_f_array[jj][2];
                }

                noface = noface + n_faces;
  
            }
            /* Populate kd tree at node zero*/
            for (jj = 0; jj<noface; jj++){
                assert(kd_insert3(kd, x_f_array_total[jj][0], x_f_array_total[jj][1], x_f_array_total[jj][2], 0) == 0);
                tree_count++;
            }
            Message0("Number of faces added to tree: %d at node %d\n" , tree_count, myid);
        }
    #endif

    /* Send array of all face coordinates to all non-zero nodes */
    #if RP_NODE
        if(I_AM_NODE_ZERO_P)
        {
            compute_node_loop_not_zero(k){
                PRF_CSEND_INT(k, &noface, 1, myid);
                PRF_CSEND_REAL(k, x_f_array_total[0], ND_ND*noface, myid);
            }
        }
    #endif

    /* Collect array of all face coordinates at all non-zero nodes */
    #if RP_NODE
        if(! I_AM_NODE_ZERO_P)
        {
            PRF_CRECV_INT(node_zero, &noface, 1, node_zero);
            x_f_array_total = (real (*)[ND_ND])malloc(ND_ND*noface*sizeof(real));
            PRF_CRECV_REAL(node_zero, x_f_array_total[0], ND_ND*noface, node_zero);
        }
    #endif


    /* Create KD tree for all non-zero compute nodes */
    #if RP_NODE
        if(! I_AM_NODE_ZERO_P)
        {
            for (jj = 0; jj<noface; jj++){
                assert(kd_insert3(kd, x_f_array_total[jj][0], x_f_array_total[jj][1], x_f_array_total[jj][2], 0) == 0);
                tree_count++;
            } 
        Message("Number of faces added to tree: %d at node %d\n" , tree_count, myid);
        }
    #endif
    
    free(x_f_array);
    free(x_f_array_total);
    PRF_GSYNC();

    /* Main Loop - Apply parameter value to each cell */

    int nocells = 0;
    begin_c_loop_int(c, vol_thread) /* Loop over all cells in each compute node */
    {
        nocells++;
        C_CENTROID(x, c, vol_thread); 
        
        /* Find closest face in kd tree and compute distance*/
        kd_result = kd_nearest3(kd, x[0], x[1], x[2]);
        face_int = kd_res_item(kd_result, pos);
        dist = sqrt(pow(x[0] - pos[0],2) + pow(x[1] - pos[1],2) + pow(x[2] - pos[2],2) );

        /* Set parameter value in cell*/
        if (dist < dist_threshold)
            phi = wall_porosity + slope * dist;
        else
            phi = core_porosity;

        perm = pow(250./2*0.000001, 2) / (8.*(1-phi)) * (-log(1.-phi) - (1-pow(1 - phi, 2))/(1+pow(1-phi, 2)) ); 

        c2 = 2. / sqrt(perm) * 3 * pow(1 - phi, 2);
        C_PROFILE(c, vol_thread, i) = c2;
        C_UDMI(c, vol_thread, 2) = c2;

        kd_res_free(kd_result);
    }
    end_c_loop_int(c, vol_thread)

    kd_free(kd);
    #if RP_NODE
        Message("Iterated over %d cells for node: %d\n", nocells, myid);
    #endif /* RP_NODE */
 }
#endif
}


