#include <string>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFab.H>
#include <hdf5.h>

void WriteMultiFabToHDF5(const std::string& plotfile, const std::string& hdf5file) {

    const amrex::PlotFileData plotfile_data(plotfile);
    const amrex::Vector<std::string> varnames = plotfile_data.varNames();
    const double time = plotfile_data.time();
    const int level = 0; // Assume there is only one level, we are reading from level 0
    const int time_step = plotfile_data.levelStep(level);
    const double dx = plotfile_data.cellSize(level)[0];
    const double dy = plotfile_data.cellSize(level)[1];

    amrex::MultiFab mf;

    amrex::VisMF::Read(mf, plotfile+"/Level_0/Cell");

    const amrex::BoxArray ba = mf.boxArray();

    // Get the dimensions of the entire MultiFab
    // Expecting a cell centered box with low and high index bounds {0,0} to {nx-1,ny-1}
    const amrex::Box minimal_box = ba.minimalBox();
    AMREX_ASSERT(minimal_box.smallEnd(0) == 0);
    AMREX_ASSERT(minimal_box.smallEnd(1) == 0);
    const int nx = minimal_box.length(0);
    const int ny = minimal_box.length(1);

    const hid_t file_id = H5Fcreate(hdf5file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    {
      const hid_t attr_space_id = H5Screate(H5S_SCALAR);
      const hid_t attr_id = H5Acreate(file_id, "time", H5T_NATIVE_DOUBLE, attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &time);
      H5Aclose(attr_id);
      H5Sclose(attr_space_id);
    }

    {
        const hid_t attr_space_id = H5Screate(H5S_SCALAR);
        const hid_t attr_id = H5Acreate(file_id, "time_step", H5T_NATIVE_INT, attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_id, H5T_NATIVE_INT, &time_step);
        H5Aclose(attr_id);
        H5Sclose(attr_space_id);
    }

    {
      const hid_t attr_space_id = H5Screate(H5S_SCALAR);
      const hid_t attr_id = H5Acreate(file_id, "dx", H5T_NATIVE_DOUBLE, attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &dx);
      H5Aclose(attr_id);
      H5Sclose(attr_space_id);
    }

    {
      const hid_t attr_space_id = H5Screate(H5S_SCALAR);
      const hid_t attr_id = H5Acreate(file_id, "dy", H5T_NATIVE_DOUBLE, attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &dy);
      H5Aclose(attr_id);
      H5Sclose(attr_space_id);
    }

    {
        const hid_t attr_space_id = H5Screate(H5S_SCALAR);
        const hid_t attr_id = H5Acreate(file_id, "nx", H5T_NATIVE_INT, attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_id, H5T_NATIVE_INT, &nx);
        H5Aclose(attr_id);
        H5Sclose(attr_space_id);
    }

    {
        const hid_t attr_space_id = H5Screate(H5S_SCALAR);
        const hid_t attr_id = H5Acreate(file_id, "ny", H5T_NATIVE_INT, attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_id, H5T_NATIVE_INT, &ny);
        H5Aclose(attr_id);
        H5Sclose(attr_space_id);
    }

    // Add an attribute to specify the data layout of the following datasets (row-major or column-major)
    {
        const char* layout = "row-major";
        const hid_t str_type = H5Tcopy(H5T_C_S1);
        H5Tset_size(str_type, strlen(layout) + 1); // +1 for the null terminator
    
        const hid_t attr_space_id = H5Screate(H5S_SCALAR);
        const hid_t attr_id = H5Acreate(file_id, "data_layout", str_type, attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr_id, str_type, layout);
        H5Aclose(attr_id);
        H5Sclose(attr_space_id);
        H5Tclose(str_type);        
    }

    // Iterate over the components of the MultiFab
    for (int component_idx = 0; component_idx < mf.nComp(); ++component_idx) {
        // Create a dataset for this component
        const hsize_t dims[2] = {static_cast<hsize_t>(nx), static_cast<hsize_t>(ny)};
        const hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
        const hid_t dataset_id = H5Dcreate(file_id, varnames[component_idx].c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        std::vector<double> component_data(dims[0]*dims[1]);

        // Loop over all sub-boxes in the MultiFab and fill component_data
        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();
            const amrex::Dim3 lo = amrex::lbound(bx);
            const amrex::Dim3 hi = amrex::ubound(bx);
            const amrex::Array4<const amrex::Real> &fab = mf.array(mfi);
            std::vector<double> data(bx.numPts());

            // Looping over the box in column-major order since that is how the multi-fab is stored
            for (int j = lo.y; j <= hi.y; ++j) {
              for (int i = lo.x; i <= hi.x; ++i) {
                    const int idx = i*dims[1] + j;  // Writing data to hdf5 file in row-major order
                    component_data[idx] = fab(i, j, 0, component_idx);
                }
            }
        }

        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, component_data.data());

        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }

    H5Fclose(file_id);
}

int main(int argc, char* argv[]) {
    amrex::Initialize(argc, argv);

    // Read in parameters from inputs file
    amrex::ParmParse pp;
    
    // Read in plotfile name
    std::string plotfile;
    pp.query("infile", plotfile);
    if (plotfile.empty()) {
        amrex::Abort("You must specify `infile'");
    }

    // Read in optional output hdf5 file name
    // Output hdf5 file (default to input plotfile name with .h5 extension)
    std::string hdf5file = plotfile + ".h5";
    pp.query("outfile", hdf5file);
    
    amrex::Print() << "Input Plotfile: " << plotfile << std::endl;

    WriteMultiFabToHDF5(plotfile, hdf5file);

    amrex::Print() << "Output HDF5 file: " << hdf5file << std::endl;

    amrex::Finalize();
    return 0;
}