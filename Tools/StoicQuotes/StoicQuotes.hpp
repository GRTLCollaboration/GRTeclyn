/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef STOIC_QUOTES_HPP_
#define STOIC_QUOTES_HPP_

#include <string_view>

namespace StoicQuotes
{
/// Select a random quote appropriate to whether the job succeeded.
[[nodiscard]] std::string_view random_quote(bool succeeded);
} // namespace StoicQuotes

#endif /* STOIC_QUOTES_HPP_ */
