#ifndef SWM_MINI_APP_KERNELS_H_
#define SWM_MINI_APP_KERNELS_H_

#include <cmath>
#include <AMReX.H>

#include "swm_mini_app_utils.h"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void UpdateIntermediateVariablesKernel( const int i, const int j, const int k,
                                        const double fsdx, const double fsdy,
                                        const amrex::Array4<amrex::Real const>& p,
                                        const amrex::Array4<amrex::Real const>& u,
                                        const amrex::Array4<amrex::Real const>& v,
                                        const amrex::Array4<amrex::Real>& cu,
                                        const amrex::Array4<amrex::Real>& cv,
                                        const amrex::Array4<amrex::Real>& h,
                                        const amrex::Array4<amrex::Real>& z)
{
    cu(i,j,k) = 0.5*(p(i,j,k) + p(i+1,j,k))*u(i,j,k);
    cv(i,j,k) = 0.5*(p(i,j,k) + p(i,j+1,k))*v(i,j,k);
    z(i,j,k) = (fsdx*(v(i+1,j,k)-v(i,j,k)) - fsdy*(u(i,j+1,k)-u(i,j,k)))/(p(i,j,k)+p(i+1,j,k)+p(i,j+1,k)+p(i+1,j+1,k));
    h(i,j,k) = p(i,j,k) + 0.25*(u(i-1,j,k)*u(i-1,j,k) + u(i,j,k)*u(i,j,k) + v(i,j-1,k)*v(i,j-1,k) + v(i,j,k)*v(i,j,k));
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void UpdateNewVariablesKernel( const int i, const int j, const int k, 
                               const double tdtsdx, const double tdtsdy, const double tdts8,
                               const amrex::Array4<amrex::Real const>& p_old,
                               const amrex::Array4<amrex::Real const>& u_old,
                               const amrex::Array4<amrex::Real const>& v_old,
                               const amrex::Array4<amrex::Real const>& cu,
                               const amrex::Array4<amrex::Real const>& cv,
                               const amrex::Array4<amrex::Real const>& h,
                               const amrex::Array4<amrex::Real const>& z,
                               const amrex::Array4<amrex::Real>& p_new,
                               const amrex::Array4<amrex::Real>& u_new,
                               const amrex::Array4<amrex::Real>& v_new)
{
    u_new(i,j,k) = u_old(i,j,k) + tdts8 * (z(i,j-1,k)+z(i,j,k)) * (cv(i,j-1,k) + cv(i,j,k) + cv(i+1,j-1,k) + cv(i+1,j,k)) - tdtsdx * (h(i+1,j,k) - h(i,j,k));
    v_new(i,j,k) = v_old(i,j,k) - tdts8 * (z(i-1,j,k)+z(i,j,k)) * (cu(i-1,j,k) + cu(i-1,j+1,k) + cu(i,j,k) + cu(i,j+1,k)) - tdtsdy * (h(i,j+1,k) - h(i,j,k));
    p_new(i,j,k) = p_old(i,j,k) - tdtsdx * (cu(i,j,k) - cu(i-1,j,k)) - tdtsdy * (cv(i,j,k) - cv(i,j-1,k));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void UpdateOldVariablesKernel( const int i, const int j, const int k, 
                               const double alpha,
                               const amrex::Array4<amrex::Real const>& p,
                               const amrex::Array4<amrex::Real const>& u,
                               const amrex::Array4<amrex::Real const>& v,
                               const amrex::Array4<amrex::Real const>& p_new,
                               const amrex::Array4<amrex::Real const>& u_new,
                               const amrex::Array4<amrex::Real const>& v_new,
                               const amrex::Array4<amrex::Real>& p_old,
                               const amrex::Array4<amrex::Real>& u_old,
                               const amrex::Array4<amrex::Real>& v_old)
{
    amrex::Real u_old_temp = u_old(i,j,k);
    amrex::Real v_old_temp = v_old(i,j,k);
    amrex::Real p_old_temp = p_old(i,j,k);

    u_old(i,j,k) = u(i,j,k) + alpha * (u_new(i,j,k) - 2.0*u(i,j,k) + u_old_temp);
    v_old(i,j,k) = v(i,j,k) + alpha * (v_new(i,j,k) - 2.0*v(i,j,k) + v_old_temp);
    p_old(i,j,k) = p(i,j,k) + alpha * (p_new(i,j,k) - 2.0*p(i,j,k) + p_old_temp);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void UpdateNewVariablesKernel( const int i, const int j, const int k, 
                               const double dx, const double dy, const double dt,
                               const amrex::Array4<amrex::Real const>& p_old,
                               const amrex::Array4<amrex::Real const>& u_old,
                               const amrex::Array4<amrex::Real const>& v_old,
                               const amrex::Array4<amrex::Real>& p_new,
                               const amrex::Array4<amrex::Real>& u_new,
                               const amrex::Array4<amrex::Real>& v_new)
{
    amrex::Real du_dt = 0.0;
    //{
    //    // u_rhs is evaluated at the y face i,j (same location as u)... 
    //    // So:
    //    //   terms with d()/dx need the difference inside of () to be evaluated at the nodes (i+1,j) and (i,j)
    //    //   terms with d()/dy need the difference inside of () to be evaluated at the cell centers (i,j) and (i,j-1)

    //    const amrex::Real dv_dx = ( InterpolateXFaceToNode(v_old,i+1,j,k) - InterpolateXFaceToNode(v_old,i,j,k)) / dx;

    //    const amrex::Real du_dy = ( InterpolateYFaceToCellCenter(u_old,i,j,k)   - InterpolateYFaceToCellCenter(u_old,i,j-1,k)) / dy;

    //    const amrex::Real v = InterpolateXFaceToYFace(v_old,i,j,k);

    //    // The pressure plus kinetic energy term is evaluated at nodes (left and right of the y face value at i,j,k)
    //    const amrex::Real PKe_r = p_old(i+1,j,k) + 0.5*(pow(InterpolateXFaceToNode(v_old,i+1,j,k),2.0) + pow(InterpolateYFaceToNode(u_old,i+1,j,k),2.0) );
    //    const amrex::Real PKe_l = p_old(i,j,k)   + 0.5*(pow(InterpolateXFaceToNode(v_old,i,j,k),2.0)   + pow(InterpolateYFaceToNode(u_old,i,j,k),  2.0) );
    //    // The d()/dx of the pressure plus kinetic energy term is evaluated at y face 
    //    dPKe_dx = (PKe_r - PKe_l) / dx;

    //    du_dt += (dv_dx - du_dy)*v - dPKe_dx;
    //}
    {
        // u_rhs is evaluated at the y face i,j (same location as u)... 
        // So all indexing in this block is relative to the y face i,j:
        //   terms with d()/dx need the difference inside of () to be evaluated at the nodes (i+1,j) and (i,j)
        //   terms with d()/dy need the difference inside of () to be evaluated at the cell centers (i,j) and (i,j-1)

        const amrex::Real u = u_old(i,j,k);
        const amrex::Real u_node_left   = InterpolateYFaceToNode(u_old,i,j,k);
        const amrex::Real u_node_right  = InterpolateYFaceToNode(u_old,i+1,j,k);
        const amrex::Real u_cell_top    = InterpolateYFaceToCellCenter(u_old,i,j,k);
        const amrex::Real u_cell_bottom = InterpolateYFaceToCellCenter(u_old,i,j-1,k);

        const amrex::Real v             = InterpolateYFaceToXFace(v_old,i,j,k);
        const amrex::Real v_node_left   = InterpolateYFaceToNode(v_old,i,j,k);
        const amrex::Real v_node_right  = InterpolateYFaceToNode(v_old,i+1,j,k);
        const amrex::Real v_cell_top    = InterpolateYFaceToCellCenter(v_old,i,j,k);
        const amrex::Real v_cell_bottom = InterpolateYFaceToCellCenter(v_old,i,j-1,k);

        const amrex::Real p             = InterpolateNodeToXFace(p_old,i,j,k);
        const amrex::Real p_node_left   = p_old(i,j,k);
        const amrex::Real p_node_right  = p_old(i+1,j,k);
        const amrex::Real p_cell_top    = InterpolateNodeToCellCenter(p_old,i,j,k);
        const amrex::Real p_cell_bottom = InterpolateNodeToCellCenter(p_old,i,j-1,k);

        const amrex::Real dv_dx = (v_node_right - v_node_left) / dx;

        const amrex::Real du_dy = (u_cell_top - u_cell_bottom) / dy;

        // The pressure plus kinetic energy term is evaluated at nodes (left and right of the y face value at i,j,k)
        const amrex::Real PKe_node_right = p_node_right + 0.5*( pow(u_node_right,2.0) + pow(v_node_right,2.0) );
        const amrex::Real PKe_node_left  = p_node_left  + 0.5*( pow(u_node_left ,2.0) + pow(v_node_left ,2.0) );
        const amrex::Real dPKe_dx = (PKe_node_right - PKe_node_left) / dx;

        du_dt += (dv_dx - du_dy)*v - dPKe_dx;
    }

    amrex::Real dv_dt = 0.0;
    {
        // v_rhs is evaluated at the x face i,j (same location as v)... 
        // So all indexing in this block is relative to the v face i,j:
        //   terms with d()/dx need the difference inside of () to be evaluated at the cell centers (i-1,j) and (i,j)
        //   terms with d()/dy need the difference inside of () to be evaluated at the nodes (i,j+1) and (i,j)

        const amrex::Real u = InterpolateYFaceToXFace(u_old,i,j,k);
        const amrex::Real u_node_top     = InterpolateYFaceToNode(u_old,i,j+1,k);
        const amrex::Real u_node_bottom  = InterpolateYFaceToNode(u_old,i,j,k);
        const amrex::Real u_cell_right   = InterpolateYFaceToCellCenter(u_old,i,j,k);
        const amrex::Real u_cell_left    = InterpolateYFaceToCellCenter(u_old,i-1,j,k);

        const amrex::Real v             = v_old(i,j,k);
        const amrex::Real v_node_top    = InterpolateYFaceToNode(v_old,i,j+1,k);
        const amrex::Real v_node_bottom = InterpolateYFaceToNode(v_old,i,j,k);
        const amrex::Real v_cell_right  = InterpolateYFaceToCellCenter(v_old,i,j,k);
        const amrex::Real v_cell_left   = InterpolateYFaceToCellCenter(v_old,i-1,j,k);

        const amrex::Real p             = InterpolateNodeToYFace(p_old,i,j,k);
        const amrex::Real p_node_top    = p_old(i,j+1,k);
        const amrex::Real p_node_bottom = p_old(i,j,k);
        const amrex::Real p_cell_right  = InterpolateNodeToCellCenter(p_old,i,j,k);
        const amrex::Real p_cell_left   = InterpolateNodeToCellCenter(p_old,i-1,j,k);

        const amrex::Real dv_dx = (v_cell_right - v_cell_left) / dx;

        const amrex::Real du_dy = (u_node_top - u_node_bottom) / dy;

        // The pressure plus kinetic energy term is evaluated at nodes (left and right of the y face value at i,j,k)
        const amrex::Real PKe_node_top    = p_node_top     + 0.5*( pow(u_node_top    ,2.0) + pow(v_node_top    ,2.0) );
        const amrex::Real PKe_node_bottom = p_node_bottom  + 0.5*( pow(u_node_bottom ,2.0) + pow(v_node_bottom ,2.0) );
        const amrex::Real dPKe_dy = (PKe_node_top - PKe_node_bottom) / dy;

        dv_dt += (dv_dx - du_dy)*u - dPKe_dy;
    }

    amrex::Real dp_dt = 0.0;
    {
        // p_rhs is evaluated at the node i,j (same location as p)... 
        // So all indexing in this block is relative to the node i,j:
        //   terms with d()/dx need the difference inside of () to be evaluated at the y faces (i-1,j) and (i,j)
        //   terms with d()/dy need the difference inside of () to be evaluated at the x faces (i,j-1) and (i,j)

        const amrex::Real u_y_face_right = u_old(i,j,k);
        const amrex::Real u_y_face_left  = u_old(i-1,j,k);

        const amrex::Real v_x_face_top    = v_old(i,j,k);
        const amrex::Real v_x_face_bottom = v_old(i,j-1,k);

        const amrex::Real p_x_face_top    = InterpolateNodeToXFace(p_old,i,j,k);
        const amrex::Real p_x_face_bottom = InterpolateNodeToXFace(p_old,i,j-1,k);
        const amrex::Real p_y_face_right  = InterpolateNodeToYFace(p_old,i,j,k);
        const amrex::Real p_y_face_left   = InterpolateNodeToYFace(p_old,i-1,j,k);

        // The pressure plus kinetic energy term is evaluated at nodes (left and right of the y face value at i,j,k)
        const amrex::Real pu_y_face_right = p_y_face_right*u_y_face_right;
        const amrex::Real pu_y_face_left  = p_y_face_left*u_y_face_left;
        const amrex::Real dpu_dx = (pu_y_face_right - pu_y_face_left) / dx;

        const amrex::Real pv_x_face_top    = p_x_face_top*v_x_face_top;
        const amrex::Real pv_x_face_bottom = p_x_face_bottom*v_x_face_bottom;
        const amrex::Real dpv_dy = (pv_x_face_top - pv_x_face_bottom) / dy;

        dp_dt += (-1.0)*(dpu_dx + dpv_dy);
    }

    u_new(i,j,k) = u_old(i,j,k) + dt * du_dt;
    v_new(i,j,k) = v_old(i,j,k) + dt * dv_dt;
    p_new(i,j,k) = p_old(i,j,k) + dt * dp_dt;
}

#endif // SWM_MINI_APP_KERNELS_H_
