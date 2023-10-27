#include "TestsArgs.hpp"

namespace Tests
{
// I know this is a horrible global non-const variable but I couldn't think
// of a better way to pass this information to Catch2 test cases.
Args g_args;
} // namespace Tests
