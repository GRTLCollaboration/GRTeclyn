/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */
#ifndef INTERVAL_H_
#define INTERVAL_H_

#include <AMReX_GpuQualifiers.H>

struct Interval
{
    Interval() = default;

    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    Interval(int a_firstComp, int a_lastComp)
        : m_begin(a_firstComp), m_end(a_lastComp)
    {
    }

    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    void define(int a_firstComp, int a_lastComp)
    {
        m_begin = a_firstComp;
        m_end   = a_lastComp;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int begin() const
    {
        return m_begin;
    }

    //! return last component number
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int end() const
    {
        return m_end;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int size() const
    {
        return m_end - m_begin + 1;
    }

    [[nodiscard]] bool contains(int a_val) const
    {
        return a_val >= m_begin && a_val <= m_end;
    }

    bool operator==(const Interval &a_interval) const
    {
        return ((m_begin == a_interval.m_begin) && (m_end == a_interval.m_end));
    }

  private:
    int m_begin{0}, m_end{-1};
};

#endif
