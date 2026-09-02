/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef STOIC_SIGNAL_HANDLING_HPP_
#define STOIC_SIGNAL_HANDLING_HPP_

#include "StoicQuotes.hpp"
#include "AMReX_BLBackTrace.H"

#include <array>
#include <csignal>
#include <string>
#include <unistd.h>

namespace StoicSignalHandling
{
[[nodiscard]] inline std::string boxed_quote(bool succeeded)
{
    const std::string message =
        " Stoic wisdom: " + std::string(StoicQuotes::random_quote(succeeded)) +
        " ";
    const std::string border(message.size(), '-');

    std::string box;
    box.reserve(3 * message.size() + 9);
    box += '+';
    box += border;
    box += "+\n|";
    box += message;
    box += "|\n+";
    box += border;
    box += "+\n";
    return box;
}

[[nodiscard]] inline std::string &signal_failure_quote()
{
    static std::string quote;
    return quote;
}

inline void handler(int signal_number)
{
    const std::string &quote = signal_failure_quote();
    static_cast<void>(::write(STDERR_FILENO, quote.data(), quote.size()));
    amrex::BLBackTrace::handler(signal_number);
}

inline void install()
{
    // Select and format the quote before a signal occurs so the handler only
    // needs to write already-prepared data.
    signal_failure_quote() = boxed_quote(false);

    constexpr std::array handled_signals = {SIGSEGV, SIGTERM, SIGINT,
                                            SIGABRT, SIGFPE,  SIGILL};
    using SignalHandler                  = void (*)(int);
    const SignalHandler error_handler =
        SIG_ERR; // NOLINT(performance-no-int-to-ptr)

    for (const int signal_number : handled_signals)
    {
        const SignalHandler previous_handler =
            std::signal(signal_number, handler);
        if (previous_handler == error_handler)
        {
            continue;
        }

        if (previous_handler != amrex::BLBackTrace::handler)
        {
            // AMReX was not handling this signal, so leave it unchanged.
            std::signal(signal_number, previous_handler);
        }
    }
}
} // namespace StoicSignalHandling

#endif /* STOIC_SIGNAL_HANDLING_HPP_ */
