/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef FILESYSTEMTOOLS_HPP_
#define FILESYSTEMTOOLS_HPP_

#include <AMReX_Utility.H>

// Some filesystem useful functions

namespace FilesystemTools
{

static bool directory_exists(const std::string &path)
{
    return amrex::FileSystem::Exists(path);
}

static void mkdir_recursive(const std::string &path)
{
    amrex::UtilCreateDirectoryDestructive(path);
}
} // namespace FilesystemTools

#endif /* FILESYSTEMTOOLS_HPP_ */
