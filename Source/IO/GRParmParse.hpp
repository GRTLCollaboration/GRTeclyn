/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef GRPARMPARSE_HPP_
#define GRPARMPARSE_HPP_

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <string>
#include <string_view>
#include <vector>

class GRParmParse : public amrex::ParmParse
{
  public:
    using amrex::ParmParse::ParmParse; // Just use ParmParse's constructor

    [[nodiscard]] static bool warnings_issued() { return s_warning_issued; }

    void error(const char *name, const std::string_view a_error_message) const
    {
        const std::string error_message =
            diagnostic_message("Error", name, a_error_message);
        amrex::Abort(error_message);
    }

    void warning(const char *name,
                 const std::string_view a_warning_message) const
    {
        s_warning_issued = true;
        const std::string warning_message =
            diagnostic_message("Warning", name, a_warning_message);
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            amrex::Warning(warning_message.c_str());
        }
    }

  protected:
    inline static bool s_warning_issued{false};

    [[nodiscard]] std::string
    diagnostic_message(const char *a_diagnostic_type, const char *name,
                       const std::string_view a_message) const
    {
        std::vector<std::string> values;
        getarr(name, values);

        std::string value_string;
        if (values.size() > 1)
        {
            value_string += "[";
        }
        for (std::size_t ivalue = 0; ivalue < values.size(); ++ivalue)
        {
            if (ivalue > 0)
            {
                value_string += " ";
            }
            value_string += values[ivalue];
        }
        if (values.size() > 1)
        {
            value_string += "]";
        }

        return std::string(a_diagnostic_type) + " from parameter " +
               this->prefixedName(name) + " = " + value_string + ": " +
               std::string(a_message);
    }
};

#endif /* GRPARMPARSE_HPP_ */
