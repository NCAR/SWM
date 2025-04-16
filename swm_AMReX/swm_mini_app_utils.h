#ifndef SWM_MINI_APP_UTILS_H_
#define SWM_MINI_APP_UTILS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

void ParseInput(int & nx, int & ny,
                amrex::Real & dx, amrex::Real & dy,
                int & max_chunk_size,
                int & n_time_steps, amrex::Real & dt,
                int & plot_interval);

void InitializeGeometry(const int nx, const int ny,
                        const amrex::Real dx, const amrex::Real dy,
                        amrex::Geometry & geom);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real LinearMapCoordinates(const amrex::Real x, 
                                 const amrex::Real x_min,  const amrex::Real x_max,
                                 const amrex::Real xi_min, const amrex::Real xi_max);

void DefineCellCenteredMultiFab(const int nx, const int ny,
                                const int max_chunk_size,
                                amrex::MultiFab & cell_centered_MultiFab);

void DefineXFaceMultiFab(const amrex::MultiFab & cell_centered_MultiFab,
                         amrex::MultiFab & x_face_MultiFab);

void DefineYFaceMultiFab(const amrex::MultiFab & cell_centered_MultiFab,
                         amrex::MultiFab & y_face_MultiFab);

void DefineNodalMultiFab(const amrex::MultiFab & cell_centered_MultiFab,
                         amrex::MultiFab & nodal_MultiFab);

void InitializeVariables(const amrex::Geometry & geom,
                         amrex::MultiFab & psi,
                         amrex::MultiFab & p,
                         amrex::MultiFab & u,
                         amrex::MultiFab & v);

void WriteOutput(const amrex::MultiFab & psi,
                 const amrex::MultiFab & p,
                 const amrex::MultiFab & u,
                 const amrex::MultiFab & v,
                 const amrex::Geometry & geom,
                 const amrex::Real time,
                 const int time_step,
                 amrex::MultiFab & output_values);

amrex::MultiFab CreateMultiFab(const amrex::MultiFab & mf);

void Copy(const amrex::MultiFab & src, amrex::MultiFab & dest);

void Swap(amrex::MultiFab & src, amrex::MultiFab & dest);

void UpdateIntermediateVariables(amrex::Real dx, amrex::Real dy, const amrex::Geometry& geom,
                                 const amrex::MultiFab& p, const amrex::MultiFab& u, const amrex::MultiFab& v,
                                 amrex::MultiFab& cu, amrex::MultiFab& cv, amrex::MultiFab& h, amrex::MultiFab& z);

void UpdateNewVariables(const double dx, const double dy, const double tdt, const amrex::Geometry& geom,
                        const amrex::MultiFab& p_old, const amrex::MultiFab& u_old, const amrex::MultiFab& v_old,
                        const amrex::MultiFab& cu, const amrex::MultiFab& cv, const amrex::MultiFab& h, const amrex::MultiFab& z,
                        amrex::MultiFab& p_new, amrex::MultiFab& u_new, amrex::MultiFab& v_new);

void UpdateOldVariables(const double alpha, const int time_step, const amrex::Geometry& geom, 
                        const amrex::MultiFab& p, const amrex::MultiFab& u, const amrex::MultiFab& v, 
                        const amrex::MultiFab& p_new, const amrex::MultiFab& u_new, const amrex::MultiFab& v_new, 
                        amrex::MultiFab& p_old, amrex::MultiFab& u_old, amrex::MultiFab& v_old);

void UpdateVariables(const amrex::Geometry& geom,
                     const amrex::MultiFab& u_new, const amrex::MultiFab& v_new, const amrex::MultiFab& p_new,
                     amrex::MultiFab& u, amrex::MultiFab& v, amrex::MultiFab& p);

void UpdateVariables(const amrex::Geometry& geom,
                     amrex::MultiFab& u_new, amrex::MultiFab& v_new, amrex::MultiFab& p_new,
                     amrex::MultiFab& u, amrex::MultiFab& v, amrex::MultiFab& p);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateXFaceToNode(const amrex::Array4<amrex::Real const>& phi_x_face, 
                                   const int i, const int j, const int k);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateXFaceToYFace(const amrex::Array4<amrex::Real const>& phi_x_face, 
                                    const int i, const int j, const int k);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateXFaceToCellCenter(const amrex::Array4<amrex::Real const>& phi_x_face, 
                                         const int i, const int j, const int k);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateYFaceToNode(const amrex::Array4<amrex::Real const>& phi_y_face, 
                                   const int i, const int j, const int k);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateYFaceToXFace(const amrex::Array4<amrex::Real const>& phi_y_face, 
                                    const int i, const int j, const int k);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateYFaceToCellCenter(const amrex::Array4<amrex::Real const>& phi_y_face, 
                                         const int i, const int j, const int k);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateNodeToXFace(const amrex::Array4<amrex::Real const>& phi_node, 
                                   const int i, const int j, const int k);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateNodeToYFace(const amrex::Array4<amrex::Real const>& phi_node, 
                                   const int i, const int j, const int k);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateNodeToCellCenter(const amrex::Array4<amrex::Real const>& phi_node, 
                                        const int i, const int j, const int k);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateCellCenterToNode(const amrex::Array4<amrex::Real const>& phi_cell_center, 
                                        const int i, const int j, const int k);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateCellCenterToXFace(const amrex::Array4<amrex::Real const>& phi_cell_center, 
                                        const int i, const int j, const int k);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real InterpolateCellCenterToYFace(const amrex::Array4<amrex::Real const>& phi_cell_center, 
                                        const int i, const int j, const int k);

void UpdateNewVariables(const double dx, const double dy, const double dt,
                        const amrex::MultiFab& p_old, const amrex::MultiFab& u_old, const amrex::MultiFab& v_old,
                        amrex::MultiFab& p_new, amrex::MultiFab& u_new, amrex::MultiFab& v_new);

#endif // SWM_MINI_APP_UTILS_H_
