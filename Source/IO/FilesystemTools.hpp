/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef FILESYSTEMTOOLS_HPP_
#define FILESYSTEMTOOLS_HPP_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

// Some filesystem useful functions

namespace FilesystemTools
{

// This is a collective operation: every MPI rank must call it
static void ensure_directory_exists(const std::string &path)
{
    if (!amrex::FileSystem::Exists(path))
    {
        // The existence check need not agree between ranks because the
        // barrier in this conditional call is disabled
        amrex::UtilCreateDirectoryDestructive(path, false);
    }

    // Do not let other ranks use the directory before it has been created
    amrex::ParallelDescriptor::Barrier(
        "FilesystemTools::ensure_directory_exists");
}
} // namespace FilesystemTools

#endif /* FILESYSTEMTOOLS_HPP_ */
