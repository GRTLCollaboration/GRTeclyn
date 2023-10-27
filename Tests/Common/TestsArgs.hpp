#ifndef COMMONARGS_HPP_
#define COMMONARGS_HPP_

namespace Tests
{
class Args
{
    int m_argc;
    char **m_argv;

  public:
    Args() {}

    void set(int argc, char **argv)
    {
        m_argc = argc;
        m_argv = argv;
    }

    int get_argc() { return m_argc; }

    char **get_argv() { return m_argv; }
};

// I know this is a horrible global non-const variable but I couldn't think
// of a better way to pass this information to Catch2 test cases.
extern Args g_args;
} // namespace Tests

#endif /* COMMONARGS_HPP_ */