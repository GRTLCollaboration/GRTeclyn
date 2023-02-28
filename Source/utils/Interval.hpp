#ifndef INTERVAL_H_
#define INTERVAL_H_

struct Interval
{
    Interval() : m_begin(0), m_end(-1) {}

    Interval(int a_firstComp, int a_lastComp)
        : m_begin(a_firstComp), m_end(a_lastComp)
    {
    }

    void define(int a_firstComp, int a_lastComp)
    {
        m_begin = a_firstComp;
        m_end   = a_lastComp;
    }

    int begin() const { return m_begin; }

    //! return last component number
    int end() const { return m_end; }

    int size() const { return m_end - m_begin + 1; }

    bool contains(int a_val) const
    {
        return a_val >= m_begin && a_val <= m_end;
    }

    bool operator==(const Interval &p) const
    {
        return ((m_begin == p.m_begin) && (m_end == p.m_end));
    }

  private:
    int m_begin, m_end;
};

#endif
